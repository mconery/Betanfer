'''
count_variant_pairs_and_null_traits.py

The purpose of this script is to identify the set of variant pairs for which
ld will need to be calculated and figure out which traits can be used to 
calculate each. 

***Run in Conda finemap environment***
'''

import zarr
import numcodecs
from multiprocessing import Pool, cpu_count, Lock
from collections import defaultdict
import csv
import sys
import getopt
import itertools

# Configure number of processes
N_PROCESSES = cpu_count()

def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
--loci => Input Loci CSV file REQUIRED
--variants => Matrix of variants and traits for which they are in sig loci REQUIRED
--output => Output file REQUIRED
--max_length => Maximum length of locus to consider (in bp; default: 5Mb) OPTIONAL
-h => help""")
    sys.exit(exit_num)

###############################################################################
#############################  SET TESTING LOCATIONS ##########################
###############################################################################

loci_file="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.master.csv"
variants_file="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.trait_variant_identities.csv"
max_length=1e8
output_file="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.unique_pairs.zarr"

###############################################################################
################## CREATE PARALLELIZED PROCESSING FUNCTIONS  ##################
###############################################################################

#Create function to process reading in of a chunk of the variant/trait matrix
def process_chunk(args):
    lines_chunk, traits, loci_map = args
    local_locus_variants = defaultdict(set)
    for line in lines_chunk:
        parts = line.strip().split(',')
        var_id = parts[0]
        chr, pos_str = var_id.split(':')
        pos = int(pos_str)
        
        for i, val in enumerate(parts[1:]):
            if val == '1':
                trait = traits[i]
                key = (trait, chr)
                
                for start, end in loci_map.get(key, []):
                    if start <= pos <= end:
                        locus_key = (trait, chr, start, end)
                        local_locus_variants[locus_key].add(var_id)
    return local_locus_variants

###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main():
    # Initialize variables
    loci_file = variants_file = max_length = output_file = None
    #Collect input parameters
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "l:v:m:o:",
            ["loci=", "variants=", "max_length=", "output="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-l", "--loci"):
            loci_file = arg
        elif opt in ("-v", "--variants"):
            variants_file = arg
        elif opt in ("-m", "--max_length"):
            max_length = int(arg)
        elif opt in ("-o", "--output"):
            output_file = arg

    # Validate required arguments
    required = [loci_file, variants_file, max_length, output_file]
    if any(r is None for r in required):
        print("ERROR: Missing required arguments")
        sys.exit(2)

    # Read in list of loci and identify loci under max length threshold
    loci_map = defaultdict(list)
    with open(loci_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            chr = row['CHR']
            start = int(row['START'])
            end = int(row['END'])
            trait = row['TRAIT']
            length = int(row['LENGTH'])
            if length <= max_length:
                loci_map[(trait, chr)].append((start, end))
    
    # Parallel processing setup for reading in variant/trait matrix
    with open(variants_file) as f:
        header = next(f).strip().split(',')
        traits = header[1:]
        all_lines = list(f)
    
    # Configure parallel processing of variant/trait matrix
    num_cores = cpu_count()
    chunk_size = len(all_lines) // (num_cores * 4)  # Adjust chunking factor
    chunks = (all_lines[i:i+chunk_size] for i in range(0, len(all_lines), chunk_size))

    # Process chunks of variant/trait in parallel
    locus_variants = defaultdict(set)
    with Pool(processes=num_cores) as pool:
        tasks = ((chunk, traits, loci_map) for chunk in chunks)
        for result in pool.imap_unordered(process_chunk, tasks, chunksize=20):
            for key, variants in result.items():
                locus_variants[key].update(variants)

    # Initialize Zarr storage
    store = zarr.DirectoryStore(args.output)
    root = zarr.group(store=store, overwrite=True)
    pairs = root.zeros('pairs', shape=(0,), 
                      dtype=[('var1', 'U50'), ('var2', 'U50')],
                      chunks=(1000000,),
                      compressor=numcodecs.Blosc(cname='zstd', clevel=3))

    #Define functions needed for parallel processing of pair counting
    def init_pool(zarr_lock, zarr_store):
        global global_lock, global_pairs
        global_lock = zarr_lock
        root = zarr.group(store=zarr_store)
        global_pairs = root['pairs']

    def process_locus(locus_key):
        variants = sorted(locus_variants[locus_key])
        pairs = []
        n = len(variants)
        for i in range(n):
            for j in range(i+1, n):
                pairs.append((variants[i], variants[j]))
        return pairs
    
    # Parallel processing setup
    lock = Lock()
    loci_keys = list(locus_variants.keys())
    chunk_size = len(loci_keys) // (cpu_count() * 4)
    
    with Pool(initializer=init_pool, initargs=(lock, store),
             processes=cpu_count()) as pool:
        results = pool.imap_unordered(process_locus, loci_keys, chunksize=chunk_size)
        
        for batch in results:
            with lock:
                pairs.append([(p[0], p[1]) for p in batch])

    # Deduplicate pairs using Zarr's built-in capabilities
    unique_pairs = set()
    for pair in pairs:
        unique_pairs.add((pair['var1'], pair['var2']))
    
    # Write final unique pairs
    root.zeros('unique_pairs', shape=(len(unique_pairs),),
              dtype=[('var1', 'U50'), ('var2', 'U50')],
              chunks=(1000000,),
              compressor=numcodecs.Blosc(cname='zstd', clevel=3))
    root['unique_pairs'][:] = list(unique_pairs)

    with open(output_file, 'w') as f:
        f.write(f"Unique variant pairs: {pair_count}\n")

if __name__ == '__main__':
    main()