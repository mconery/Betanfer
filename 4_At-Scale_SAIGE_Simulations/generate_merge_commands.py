#!/usr/bin/env python3
"""
LD Matrix Merge Command Generator
Generates commands to merge AFR and EUR LD matrices using the combine_ld_matrices.py script.

This script scans directories for matching AFR and EUR matrix files and creates
commands to merge them into weighted combined matrices.

Directory structure expected:
    input_dir/
    ├── AFR/
    │   ├── matrix1.unphased.vcor1
    │   ├── matrix1.unphased.vcor1.vars
    │   ├── matrix2.unphased.vcor1
    │   └── matrix2.unphased.vcor1.vars
    ├── EUR/
    │   ├── matrix1.unphased.vcor1
    │   ├── matrix1.unphased.vcor1.vars
    │   ├── matrix2.unphased.vcor1
    │   └── matrix2.unphased.vcor1.vars
    └── weighted/
        └── [output files will be created here]

Usage:
    python generate_merge_commands.py \
        --input-dirs /path/to/dir1 /path/to/dir2 ... \
        --output-commands commands.txt \
        --afr-sample-size 10000 \
        --eur-sample-size 15000 \
        --combiner-script /path/to/combine_ld_matrices.py \
        --threads 32

"""

import argparse
import os
from pathlib import Path
import logging
from typing import List, Tuple, Dict
import glob

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MergeCommandGenerator:
    def __init__(self, input_dirs: List[str], output_commands: str, 
                 afr_sample_size: int, eur_sample_size: int,
                 combiner_script: str, threads: int = 32, 
                 use_gpu: bool = False, verbose: bool = False):
        """
        Initialize the merge command generator.

        Args:
            input_dirs: List of input directories containing AFR/EUR/weighted subdirs
            output_commands: Path to output file containing merge commands
            afr_sample_size: Sample size for AFR populations
            eur_sample_size: Sample size for EUR populations  
            combiner_script: Path to the combine_ld_matrices.py script
            threads: Number of threads to use for each merge operation
            use_gpu: Whether to use GPU acceleration
            verbose: Enable verbose logging
        """
        self.input_dirs = [Path(d) for d in input_dirs]
        self.output_commands = Path(output_commands)
        self.afr_sample_size = afr_sample_size
        self.eur_sample_size = eur_sample_size
        self.combiner_script = Path(combiner_script)
        self.threads = threads
        self.use_gpu = use_gpu
        self.verbose = verbose

        # Validate inputs
        self._validate_inputs()

        logger.info(f"Initialized command generator for {len(input_dirs)} directories")
        logger.info(f"AFR sample size: {afr_sample_size}")
        logger.info(f"EUR sample size: {eur_sample_size}")
        logger.info(f"Threads per job: {threads}")
        logger.info(f"GPU acceleration: {'Enabled' if use_gpu else 'Disabled'}")

    def _validate_inputs(self):
        """Validate input parameters and directories."""
        # Check combiner script exists
        if not self.combiner_script.exists():
            raise FileNotFoundError(f"Combiner script not found: {self.combiner_script}")

        # Check input directories exist and have required structure
        for input_dir in self.input_dirs:
            if not input_dir.exists():
                raise FileNotFoundError(f"Input directory not found: {input_dir}")

            # Check for AFR, EUR, and weighted subdirectories
            afr_dir = input_dir / "AFR"
            eur_dir = input_dir / "EUR"
            weighted_dir = input_dir / "weighted"

            if not afr_dir.exists():
                raise FileNotFoundError(f"AFR subdirectory not found: {afr_dir}")
            if not eur_dir.exists():
                raise FileNotFoundError(f"EUR subdirectory not found: {eur_dir}")
            if not weighted_dir.exists():
                logger.info(f"Creating weighted output directory: {weighted_dir}")
                weighted_dir.mkdir(parents=True, exist_ok=True)

        # Create output directory for commands file
        self.output_commands.parent.mkdir(parents=True, exist_ok=True)

    def find_matrix_files(self, directory: Path) -> Dict[str, Path]:
        """
        Find all LD matrix files in a directory.

        Returns:
            Dict mapping base filename to full path of .unphased.vcor1 file
        """
        matrix_files = {}

        # Find all .unphased.vcor1 files
        pattern = directory / "*.unphased.vcor1"
        for matrix_file in glob.glob(str(pattern)):
            matrix_path = Path(matrix_file)
            vars_file = Path(str(matrix_file) + ".vars")

            # Check that corresponding .vars file exists
            if vars_file.exists():
                # Use base filename (without .unphased.vcor1) as key
                base_name = matrix_path.name.replace('.unphased.vcor1', '')
                matrix_files[base_name] = matrix_path
            else:
                logger.warning(f"Missing .vars file for {matrix_file}")

        return matrix_files

    def find_matching_matrices(self, input_dir: Path) -> List[Tuple[str, Path, Path, Path]]:
        """
        Find matching AFR and EUR matrix files in an input directory.

        Returns:
            List of tuples: (base_name, afr_file, eur_file, output_prefix)
        """
        afr_dir = input_dir / "AFR"
        eur_dir = input_dir / "EUR"
        weighted_dir = input_dir / "weighted"

        # Find matrix files in each directory
        afr_files = self.find_matrix_files(afr_dir)
        eur_files = self.find_matrix_files(eur_dir)

        logger.info(f"Found {len(afr_files)} AFR matrices in {afr_dir}")
        logger.info(f"Found {len(eur_files)} EUR matrices in {eur_dir}")

        # Find matching pairs
        matching_pairs = []
        for base_name in afr_files.keys():
            if base_name in eur_files:
                afr_file = afr_files[base_name]
                eur_file = eur_files[base_name]
                output_prefix = weighted_dir / base_name

                matching_pairs.append((base_name, afr_file, eur_file, output_prefix))
            else:
                logger.warning(f"No matching EUR file found for AFR matrix: {base_name}")

        # Check for EUR files without AFR matches
        for base_name in eur_files.keys():
            if base_name not in afr_files:
                logger.warning(f"No matching AFR file found for EUR matrix: {base_name}")

        logger.info(f"Found {len(matching_pairs)} matching AFR-EUR matrix pairs")
        return matching_pairs

    def generate_merge_command(self, base_name: str, afr_file: Path, 
                             eur_file: Path, output_prefix: Path) -> str:
        """Generate a single merge command for an AFR-EUR matrix pair."""

        # Build the command
        cmd_parts = [
            f"python {self.combiner_script}",
            f"--ld-files {afr_file} {eur_file}",
            f"--sample-sizes {self.afr_sample_size} {self.eur_sample_size}",
            f"--output-prefix {output_prefix}",
            f"--threads {self.threads}"
        ]

        if self.use_gpu:
            # Use GPU version if requested
            gpu_script = str(self.combiner_script).replace('.py', '_gpu.py')
            cmd_parts[0] = f"python {gpu_script}"
            cmd_parts.append("--use-gpu")

        if self.verbose:
            cmd_parts.append("--verbose")

        command = ' '.join(cmd_parts)
        return command

    def generate_all_commands(self) -> List[str]:
        """Generate all merge commands for all input directories."""
        all_commands = []
        total_pairs = 0

        for input_dir in self.input_dirs:
            logger.info(f"Processing directory: {input_dir}")

            # Find matching matrix pairs in this directory
            matching_pairs = self.find_matching_matrices(input_dir)
            total_pairs += len(matching_pairs)

            # Generate commands for each pair
            for base_name, afr_file, eur_file, output_prefix in matching_pairs:
                command = self.generate_merge_command(base_name, afr_file, eur_file, output_prefix)
                all_commands.append(command)

                logger.info(f"Generated command for {base_name}")
                if self.verbose:
                    logger.info(f"  AFR: {afr_file}")
                    logger.info(f"  EUR: {eur_file}")
                    logger.info(f"  Output: {output_prefix}")

        logger.info(f"Generated {len(all_commands)} merge commands for {total_pairs} matrix pairs")
        return all_commands

    def write_commands_file(self, commands: List[str]):
        """Write all commands to the output file."""

        with open(self.output_commands, 'w') as f:
            # Write header comment
            f.write(f"#!/bin/bash\n")
            f.write(f"# Generated merge commands for AFR-EUR LD matrix combination\n")
            f.write(f"# Total commands: {len(commands)}\n")
            f.write(f"# AFR sample size: {self.afr_sample_size}\n")
            f.write(f"# EUR sample size: {self.eur_sample_size}\n")
            f.write(f"# Threads per job: {self.threads}\n")
            f.write(f"# Generated on: $(date)\n\n")

            # Write each command
            for i, command in enumerate(commands, 1):
                f.write(f"# Command {i}\n")
                f.write(f"{command}\n\n")

        # Make the file executable
        os.chmod(self.output_commands, 0o755)

        logger.info(f"Commands written to: {self.output_commands}")
        logger.info(f"File is executable - you can run: {self.output_commands}")

    def write_slurm_array_job(self, commands: List[str]):
        """Write a SLURM array job script for running all commands."""

        # Also create a simple commands-only file (no comments) for the array job
        commands_only_file = self.output_commands.parent / f"{self.output_commands.stem}.txt"
        with open(commands_only_file, 'w') as f:
            for command in commands:
                f.write(f"{command}\n")

        logger.info(f"Commands-only file for array job: {commands_only_file}")
        return commands_only_file

    def run(self):
        """Execute the complete command generation workflow."""
        logger.info("Starting merge command generation workflow...")

        try:
            # Generate all commands
            commands = self.generate_all_commands()

            if not commands:
                logger.warning("No commands generated - check input directories and file matching")
                return

            # Write commands file
            self.write_commands_file(commands)

            # Write SLURM array job script
            array_script, commands_only_file = self.write_slurm_array_job(commands)

            # Print summary
            logger.info("\n=== SUMMARY ===")
            logger.info(f"Input directories: {len(self.input_dirs)}")
            logger.info(f"Total commands generated: {len(commands)}")
            logger.info(f"Commands file: {self.output_commands}")
            logger.info(f"Array job script: {array_script}")
            logger.info(f"Commands-only file: {commands_only_file}")

            logger.info("\n=== USAGE OPTIONS ===")
            logger.info(f"1. Run all commands sequentially:")
            logger.info(f"   {self.output_commands}")
            logger.info(f"2. Submit as SLURM array job:")
            logger.info(f"   sbatch {array_script}")
            logger.info(f"3. Run individual commands manually from:")
            logger.info(f"   {self.output_commands}")

            logger.info("\nWorkflow completed successfully!")

        except Exception as e:
            logger.error(f"Workflow failed: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(
        description='Generate merge commands for AFR-EUR LD matrix pairs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-dirs', nargs='+', required=True,
                       help='List of input directories containing AFR/EUR/weighted subdirs')
    parser.add_argument('--output-commands', required=True,
                       help='Output file for generated commands')
    parser.add_argument('--afr-sample-size', type=int, required=True,
                       help='Sample size for AFR populations')
    parser.add_argument('--eur-sample-size', type=int, required=True,
                       help='Sample size for EUR populations')
    parser.add_argument('--combiner-script', required=True,
                       help='Path to combine_ld_matrices.py script')
    parser.add_argument('--threads', type=int, default=32,
                       help='Number of threads per merge operation')
    parser.add_argument('--use-gpu', action='store_true',
                       help='Use GPU-accelerated version of combiner')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate that input directories exist
    for input_dir in args.input_dirs:
        if not os.path.exists(input_dir):
            parser.error(f"Input directory not found: {input_dir}")

    # Initialize and run generator
    generator = MergeCommandGenerator(
        input_dirs=args.input_dirs,
        output_commands=args.output_commands,
        afr_sample_size=args.afr_sample_size,
        eur_sample_size=args.eur_sample_size,
        combiner_script=args.combiner_script,
        threads=args.threads,
        use_gpu=args.use_gpu,
        verbose=args.verbose
    )

    generator.run()

if __name__ == '__main__':
    main()
