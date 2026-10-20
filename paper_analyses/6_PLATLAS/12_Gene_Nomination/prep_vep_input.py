#!/usr/bin/env python3
'''
prep_vep_input.py

Builds the GRCh37 VCF that VEP annotates, from the unique fine-mapped variants
collated in step 11 (annotation is per-variant, independent of trait/ancestry).

Fine-mapping variant_id coordinates are GRCh38, but the VEP cache (release 110)
and the ABC-Max enhancer file used by nominate_genes_comprehensive.py are
GRCh37. Every unique variant_id is therefore lifted to GRCh37 with
allele-aware reference validation, and is accepted only if:
  - REF matches the GRCh38 reference at the source position;
  - the first base maps uniquely, and every other REF base maps uniquely and
    contiguously on the same strand (indels also need both flanking bases);
  - REF (reverse-complemented on minus-strand mappings) matches the GRCh37
    reference at the lifted position. Minus-strand indels are re-anchored on
    the preceding GRCh37 base so the record stays left-padded.
Rejected variants are not annotated; the QC table records why.

The exception is approved_allele_recoveries.tsv (this directory): 91 variants
that fail the checks above only because GRCh37 carries the original GRCh38 ALT
allele as its reference (71 X-chromosome SNVs confirmed by reciprocal
liftover, 20 variants in chain gaps confirmed by 25 bp of flanking sequence on
both haplotypes). They are written with the reviewed GRCh37 alleles, so the VEP
consequence is GRCh37-ALT relative; effect alleles and PIPs are unchanged.

The accepted VCF is then normalized with bcftools norm -c e (any REF mismatch
is an error) against the same GRCh37 reference and sorted, using the bcftools
bundled with pysam (the VEP environment only ships legacy bcftools 0.1.19,
which has no norm). The normalized VCF must keep every variant_id.

Input:
  --variants      gwas_finemapped_variants.tsv (step 11), columns variant_id
                  (CHR:POS:REF:ALT, GRCh38), median_pip, in_cs_any
  --recoveries    approved_allele_recoveries.tsv, columns variant_id,
                  original_rejection, recovery_method, target_chrom,
                  target_pos, target_ref, target_alt
  --source-fasta  GRCh38 FASTA (indexed); --target-fasta GRCh37 FASTA (indexed)

Output (in --out-dir):
  vep_input_all_hg19.vcf             accepted variants at GRCh37 coordinates,
                                     ID = original GRCh38 variant_id
  vep_input_all_hg19.normalized.vcf  the same records after bcftools norm and
                                     sort -- input to run_vep_annotation.sh
  vep_input_variant_qc.tsv.gz        one row per unique variant_id: status,
                                     reason, source/target coordinates, alleles

Usage:
  python prep_vep_input.py [OPTIONS]

***Run in the betainfer conda environment on Aurora (needs pysam, pyliftover;
not installed in the finemap / finemap_old environments)***
'''

import sys
import csv
import gzip
import getopt
import math
import tempfile
from collections import Counter
from pathlib import Path

import pysam
import pysam.bcftools
from pyliftover import LiftOver

HERE = Path(__file__).resolve().parent
VARIANTS_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets/gwas_finemapped_variants.tsv"
OUT_DIR_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination"
RECOVERIES_DEFAULT = str(HERE / "approved_allele_recoveries.tsv")
# Reference FASTAs are external resources only available in our trees.
SOURCE_FASTA_DEFAULT = "/lus/flare/projects/GeomicVar/alphagenome_dataset/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
TARGET_FASTA_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/interim/2026-09-23_annotation-repair/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa"

COMPLEMENT = str.maketrans("ACGT", "TGCA")
# PLINK numeric chromosome codes -> UCSC chain naming (23=X, 24=Y, 25=X [PAR], 26=MT)
CHROM_ALIASES = {"23": "X", "24": "Y", "25": "X", "26": "M", "MT": "M"}
QC_FIELDS = [
    "variant_id", "source_chrom", "source_pos", "source_ref", "source_alt",
    "nomination_eligible", "in_any_cs", "status", "reason", "detail",
    "source_ref_observed", "mapping_count", "mapping_strand", "target_chrom",
    "target_pos", "target_ref", "target_alt", "target_ref_observed",
    "reverse_indel_anchor_added", "recovery_method",
]


def help(exit_num=0):
    print("""Usage: python prep_vep_input.py [OPTIONS]
Optional:
  --variants      Fine-mapped variants table (step 11 output)
                   (default: {v})
  --recoveries    Approved GRCh37 allele recoveries (default: {r})
  --chain         hg38ToHg19 chain file (default: pyliftover's cached UCSC
                   hg38ToHg19.over.chain.gz, downloaded on first use)
  --source-fasta  Indexed GRCh38 FASTA (default: {s})
  --target-fasta  Indexed GRCh37 FASTA (default: {t})
  --out-dir       Output directory (default: {o})
  -h, --help      Show this help message""".format(
        v=VARIANTS_DEFAULT, r=RECOVERIES_DEFAULT, s=SOURCE_FASTA_DEFAULT,
        t=TARGET_FASTA_DEFAULT, o=OUT_DIR_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  REFERENCE AND LIFTOVER CHECKS  ####################################
###############################################################################

class AuditReject(Exception):
    """An expected data defect; the variant is recorded as rejected."""
    def __init__(self, reason, detail=""):
        self.reason, self.detail = reason, detail


def chain_chrom(chrom):
    name = str(chrom).removeprefix("chr")
    return "chr" + CHROM_ALIASES.get(name, name)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def resolve_contig(fasta, requested):
    """FASTA contig name for a chromosome given in any of the usual namings."""
    normalized = chain_chrom(requested)
    possibilities = [requested, normalized, normalized.removeprefix("chr")]
    if normalized == "chrM":
        possibilities.extend(["MT", "chrMT"])
    found = [name for name in dict.fromkeys(possibilities) if name in fasta.references]
    if len(found) != 1:
        reason = "missing_reference_contig" if not found else "ambiguous_reference_contig"
        raise AuditReject(reason, f"{requested}: {found}")
    return found[0]


def fetch_checked(fasta, chrom, pos, length, stage):
    """Return (contig, sequence) for a 1-based position, rejecting out-of-bounds requests."""
    try:
        contig = resolve_contig(fasta, chrom)
    except AuditReject as error:
        raise AuditReject(stage + "_" + error.reason, error.detail) from error
    try:
        if pos < 1 or pos + length - 1 > fasta.get_reference_length(contig):
            raise AuditReject(stage + "_reference_out_of_bounds")
        return contig, fasta.fetch(contig, pos - 1, pos - 1 + length).upper()
    except (ValueError, OSError, KeyError) as error:
        raise AuditReject(stage + "_reference_fetch_error", str(error)) from error


def unique_mapping(lifter, chrom, zero_pos, stage):
    """Lift one 0-based position; it must map to exactly one place."""
    mappings = lifter.convert_coordinate(chrom, zero_pos) or []
    if len(mappings) != 1:
        reason = "unmapped" if not mappings else "ambiguous_mapping"
        raise AuditReject(stage + "_" + reason, f"position={zero_pos + 1}; hits={len(mappings)}")
    dest, dest_pos, strand, _ = mappings[0]
    if strand not in ["+", "-"]:
        raise AuditReject(stage + "_invalid_strand", str(strand))
    return dest, int(dest_pos), strand


def check_interval(lifter, chrom, pos, ref, first, indel):
    """Check every REF base; indels also require both source flanks to map."""
    dest, first_pos, strand = first
    direction = 1 if strand == "+" else -1
    for offset in range(1, len(ref)):
        mapping = unique_mapping(lifter, chrom, pos - 1 + offset, "ref_span")
        if mapping != (dest, first_pos + direction * offset, strand):
            raise AuditReject("ref_span_noncontiguous", f"offset={offset}; mapped={mapping}")
    if indel:
        if pos == 1:
            raise AuditReject("indel_missing_upstream_source_flank")
        for offset in [-1, len(ref)]:
            mapping = unique_mapping(lifter, chrom, pos - 1 + offset, "indel_flank")
            if mapping != (dest, first_pos + direction * offset, strand):
                raise AuditReject("indel_flank_noncontiguous", f"offset={offset}; mapped={mapping}")


###############################################################################
##########  PER-VARIANT LIFTOVER  #############################################
###############################################################################

def lift_alleles(row, lifter, chrom, pos, ref, alt):
    """Lift a source-validated variant; returns (dest, dest_pos, ref, alt, strand)."""
    row["mapping_count"] = len(lifter.convert_coordinate(chain_chrom(chrom), pos - 1) or [])
    first = unique_mapping(lifter, chain_chrom(chrom), pos - 1, "anchor")
    row["mapping_strand"] = first[2]
    check_interval(lifter, chain_chrom(chrom), pos, ref, first, len(ref) != len(alt))
    dest, first_pos, strand = first
    dest_pos = first_pos + 1 if strand == "+" else first_pos - len(ref) + 2
    if strand == "-":
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    return dest, dest_pos, ref, alt, strand


def reanchor_reverse_indel(row, target, dest, dest_pos, ref, alt):
    """Reverse-complemented indels may have a right-hand anchor. A shared
    preceding base restores a valid left-padded representation without
    swapping alleles or claiming repeat-aware left normalization."""
    _, padding = fetch_checked(target, dest, dest_pos - 1, 1, "target_anchor")
    if padding not in "ACGT":
        raise AuditReject("target_anchor_ambiguous_base")
    dest_pos -= 1
    ref, alt = padding + ref, padding + alt
    _, observed = fetch_checked(target, dest, dest_pos, len(ref), "target_padded")
    if observed != ref:
        raise AuditReject("target_padded_ref_mismatch")
    row.update(target_pos=dest_pos, target_ref=ref, target_alt=alt,
               target_ref_observed=observed, reverse_indel_anchor_added=1)


def prepare_variant(variant_id, flags, source, target, lifter):
    """Return one complete QC row; expected data defects are explicit rejects."""
    parts = variant_id.split(":")
    row = dict.fromkeys(QC_FIELDS, "")
    row.update(variant_id=variant_id, nomination_eligible=int(bool(flags & 1)),
               in_any_cs=int(bool(flags & 2)), reverse_indel_anchor_added=0)
    try:
        if len(parts) != 4:
            raise AuditReject("invalid_variant_id")
        chrom, pos, ref, alt = parts
        row.update(source_chrom=chrom, source_pos=pos, source_ref=ref, source_alt=alt)
        try:
            pos = int(pos)
        except (ValueError, TypeError) as error:
            raise AuditReject("invalid_position", str(pos)) from error
        ref, alt = ref.upper(), alt.upper()
        if not ref or not alt or set(ref + alt) - set("ACGT"):
            raise AuditReject("unsupported_or_ambiguous_allele")
        if ref == alt:
            raise AuditReject("identical_ref_alt")
        _, row["source_ref_observed"] = fetch_checked(source, chrom, pos, len(ref), "source")
        if row["source_ref_observed"] != ref:
            raise AuditReject("source_ref_mismatch")
        dest, dest_pos, ref, alt, strand = lift_alleles(row, lifter, chrom, pos, ref, alt)
        row.update(target_chrom=dest, target_pos=dest_pos, target_ref=ref, target_alt=alt)
        dest, observed = fetch_checked(target, dest, dest_pos, len(ref), "target")
        row.update(target_chrom=dest, target_ref_observed=observed)
        if observed != ref:
            raise AuditReject("target_ref_mismatch")
        if strand == "-" and len(ref) != len(alt):
            reanchor_reverse_indel(row, target, dest, dest_pos, ref, alt)
        row.update(status="accepted", reason="pass")
    except AuditReject as error:
        row.update(status="rejected", reason=error.reason, detail=error.detail)
    return row


def apply_recovery(row, recovery, target):
    """Accept a rejected variant at its reviewed GRCh37 alleles, after re-checking REF."""
    if row["reason"] != recovery["original_rejection"]:
        print(f"ERROR: {row['variant_id']} was rejected for {row['reason']}, "
              f"but its recovery was approved for {recovery['original_rejection']}")
        sys.exit(1)
    chrom, pos, ref = recovery["target_chrom"], int(recovery["target_pos"]), recovery["target_ref"]
    if target.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
        print(f"ERROR: recovered GRCh37 REF mismatch for {row['variant_id']}")
        sys.exit(1)
    row.update(status="accepted", reason="verified_allele_aware_recovery", detail="",
               target_chrom=chrom, target_pos=pos, target_ref=ref, target_alt=recovery["target_alt"],
               target_ref_observed=ref, recovery_method=recovery["recovery_method"])


###############################################################################
##########  INPUT / OUTPUT  ###################################################
###############################################################################

def load_variants(path):
    """Collapse trait/population rows to variant_id -> flags (1 = median_pip >= 0.01, 2 = in a CS)."""
    variants = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                value = float(row["median_pip"])
                eligible = math.isfinite(value) and value >= 0.01
            except (ValueError, TypeError):
                eligible = False
            flags = int(eligible) | (2 if row["in_cs_any"] == "1" else 0)
            variants[row["variant_id"]] = variants.get(row["variant_id"], 0) | flags
    return variants


def load_recoveries(path):
    with open(path, newline="") as fh:
        return {r["variant_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def write_vcf(path, records, fasta, reference):
    """Write accepted records sorted by FASTA contig order, position, variant_id."""
    accepted = [row for row in records if row["status"] == "accepted"]
    ordering = {name: i for i, name in enumerate(fasta.references)}
    accepted.sort(key=lambda r: (ordering[r["target_chrom"]], int(r["target_pos"]), r["variant_id"]))
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##reference={Path(reference).resolve()}\n")
        for chrom in sorted({r["target_chrom"] for r in accepted}, key=ordering.get):
            fh.write(f"##contig=<ID={chrom},length={fasta.get_reference_length(chrom)}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in accepted:
            fh.write("\t".join(map(str, [r["target_chrom"], r["target_pos"], r["variant_id"],
                                         r["target_ref"], r["target_alt"], ".", "PASS", "."])) + "\n")
    return len(accepted)


def vcf_ids(path):
    with open(path) as fh:
        return [line.split("\t", 3)[2] for line in fh if not line.startswith("#")]


def normalize_vcf(raw_vcf, norm_vcf, reference):
    """bcftools norm -c e (REF mismatch = error) and sort; every variant_id must survive."""
    unsorted_vcf = norm_vcf.with_suffix(".unsorted.vcf")
    pysam.bcftools.norm("-f", str(reference), "-c", "e", "-Ov", "-o", str(unsorted_vcf),
                        str(raw_vcf), catch_stdout=False)
    print(f"NOTE: bcftools norm: {pysam.bcftools.norm.get_messages().strip()}")
    with tempfile.TemporaryDirectory(dir=norm_vcf.parent) as scratch:
        pysam.bcftools.sort("-T", scratch, "-Ov", "-o", str(norm_vcf), str(unsorted_vcf),
                            catch_stdout=False)
    unsorted_vcf.unlink()
    if sorted(vcf_ids(raw_vcf)) != sorted(vcf_ids(norm_vcf)):
        print(f"ERROR: normalization changed the variant set ({raw_vcf} -> {norm_vcf})")
        sys.exit(1)


def save_qc(path, records):
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=QC_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def driver(variants, recoveries, chain, source_fasta, target_fasta, out_dir):
    variant_flags = load_variants(variants)
    approved = load_recoveries(recoveries)
    print(f"NOTE: {len(variant_flags):,} unique variants; {len(approved)} approved recoveries")
    lifter = LiftOver(chain) if chain else LiftOver("hg38", "hg19")
    records = []
    with pysam.FastaFile(source_fasta) as source, pysam.FastaFile(target_fasta) as target:
        for index, (vid, flags) in enumerate(sorted(variant_flags.items()), 1):
            row = prepare_variant(vid, flags, source, target, lifter)
            if row["status"] == "rejected" and vid in approved:
                apply_recovery(row, approved[vid], target)
            records.append(row)
            if index % 100000 == 0:
                print(f"NOTE: checked {index:,}/{len(variant_flags):,}", flush=True)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        save_qc(out_dir / "vep_input_variant_qc.tsv.gz", records)
        n_vcf = write_vcf(out_dir / "vep_input_all_hg19.vcf", records, target, target_fasta)
    counts = Counter((r["status"], r["reason"]) for r in records)
    for (status, reason), n in sorted(counts.items()):
        print(f"NOTE: {status:>8}  {reason:<32} {n:>9,}")
    normalize_vcf(out_dir / "vep_input_all_hg19.vcf", out_dir / "vep_input_all_hg19.normalized.vcf", target_fasta)
    print(f"NOTE: {n_vcf:,} variants written to {out_dir / 'vep_input_all_hg19.normalized.vcf'}")


def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "variants=", "recoveries=", "chain=", "source-fasta=", "target-fasta=", "out-dir="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    args = dict(variants=VARIANTS_DEFAULT, recoveries=RECOVERIES_DEFAULT, chain="",
                source_fasta=SOURCE_FASTA_DEFAULT, target_fasta=TARGET_FASTA_DEFAULT,
                out_dir=OUT_DIR_DEFAULT)
    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        args[o.lstrip("-").replace("-", "_")] = a

    for key in ("variants", "recoveries", "source_fasta", "target_fasta"):
        if not Path(args[key]).is_file():
            print(f"ERROR: --{key.replace('_', '-')} not found: {args[key]}")
            sys.exit(1)
    for key in ("source_fasta", "target_fasta"):
        if not Path(args[key] + ".fai").is_file():
            print(f"ERROR: indexed FASTA required (missing {args[key]}.fai)")
            sys.exit(1)

    driver(**args)


if __name__ == "__main__":
    main()
