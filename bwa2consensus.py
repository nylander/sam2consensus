#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2023 nylander <johan.nylander@nrm.se>
#
# Distributed under terms of the MIT license.
#
# Last modified: 2026-04-08 12:39:52
#
# Sign: JN

"""

Map fastq sequences to reference using bwa mem, calculate consensus sequence
using samtools consensus.

Input: file_R1.fq.gz, file_R2.fq.gz, ref.fasta

Output: consensus.fasta

Samtools consensus documentation:
<https://www.htslib.org/doc/1.15/samtools-consensus.html>

"""

import argparse
import os
import subprocess
import sys
from typing import List, Tuple

def parse_args(args: List[str]) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Map fastq sequences to reference using bwa mem, calculate consensus sequence using samtools consensus."
    )
    parser.add_argument(
        "--fastq1",
        type=str,
        required=True,
        help="Fastq file with forward reads",
    )
    parser.add_argument(
        "--fastq2",
        type=str,
        required=True,
        help="Fastq file with reverse reads",
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference fasta file",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=False,
        default=".",
        help="Output directory",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        required=False,
        default="consensus",
        help="Prefix for output files",
    )
    parser.add_argument(
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use",
    )
    return parser.parse_args(args)

def run_bwa_index(ref: str) -> None:
    """Run bwa index."""
    # Run bwa index
    print("Running bwa index")
    cmd = f"bwa index {ref} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)

def run_samtools_index(bam: str, threads: int) -> None:
    """Run samtools index."""
    # Run samtools index
    print("Running samtools index")
    cmd = f"samtools index -@ {threads} {bam}"
    subprocess.run(cmd, shell=True, check=True)

def run_bwa(fastq1: str, fastq2: str, ref: str, outdir: str, prefix: str, threads: int) -> Tuple[str, str]:
    """Run bwa mem."""
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outbam = os.path.join(outdir, f"{prefix}.bam")
    # Run bwa mem
    print("Running bwa mem")
    cmd = f"bwa mem -t {threads} {ref} {fastq1} {fastq2} 2>/dev/null | samtools sort --threads {threads} -o {outbam} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)
    return outbam

def run_samtools_consensus(bam: str, outdir: str, prefix: str, threads: str) -> Tuple[str, str]:
    """Run samtools consensus."""
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outfasta = os.path.join(outdir, f"{prefix}.fasta")
    # Run samtools consensus
    print("Running samtools consensus")
    cmd = f"samtools consensus --threads {threads} {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_consensus_a(bam: str, outdir: str, prefix: str, threads: str) -> Tuple[str, str]:
    """
    Run samtools consensus with -a option (Outputs all bases, from start to end of reference,
    even when the aligned data does not extend to the ends).
    """
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outfasta = os.path.join(outdir, f"{prefix}.fasta")
    # Run samtools consensus
    print("Running samtools consensus with -a option")
    cmd = f"samtools consensus --threads {threads} -a {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_consensus_iupac(bam: str, outdir: str, prefix: str, threads: str) -> Tuple[str, str]:
    """Run samtools consensus with --ambig option."""
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outfasta = os.path.join(outdir, f"{prefix}.fasta")
    # Run samtools consensus
    print("Running samtools consensus with --ambig option")
    cmd = f"samtools consensus --threads {threads} --ambig {bam} -o {outfasta}"
    subprocess.run(cmd, shell=True, check=True)
    return outfasta

def run_samtools_stats(bam: str, outdir: str, prefix: str, threads: str) -> Tuple[str, str]:
    """Run samtools stats."""
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outstats = os.path.join(outdir, f"{prefix}.stats")
    # Run samtools stats
    print("Running samtools stats")
    cmd = f"samtools stats --threads {threads} {bam} > {outstats}"
    subprocess.run(cmd, shell=True, check=True)
    return outstats

def run_samtools_depth_a(bam: str, outdir: str, prefix: str, threads: str) -> Tuple[str, str]:
    """Run samtools depth with -a option"""
    # Create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Create output file name
    outdepth = os.path.join(outdir, f"{prefix}.depth")
    # Run samtools depth
    print("Running samtools depth with -a option")
    cmd = f"samtools depth -a --threads {threads} {bam} > {outdepth}"
    subprocess.run(cmd, shell=True, check=True)
    return outdepth

def main():
    """Main function."""
    args = parse_args(sys.argv[1:])
    print(f"\nStart of workflow for {args.prefix}")
    # Run run_bwa_index unless file exists with suffix .amb
    refi=args.ref + '.amb'
    if not os.path.isfile(refi):
        run_bwa_index(args.ref)
    # Run bwa index
    #run_bwa_index(args.ref)
    # Run bwa mem
    bam = run_bwa(args.fastq1, args.fastq2, args.ref, args.outdir, args.prefix, args.threads)
    # Run samtools index
    run_samtools_index(bam, args.threads)
    # Run samtools consensus
    #fasta = run_samtools_consensus(bam, args.outdir, args.prefix, args.threads)
    # Run samtools consensus with -a option
    fasta = run_samtools_consensus_a(bam, args.outdir, args.prefix, args.threads)
    # Run samtools consensus with --ambig option
    #fasta = run_samtools_consensus_iupac(bam, args.outdir, args.prefix, args.threads)
    # Run samtools depth
    #depth = run_samtools_depth(bam, args.outdir, args.prefix, args.threads)
    # Run samtools depth with -a
    depth = run_samtools_depth_a(bam, args.outdir, args.prefix, args.threads)
    # Run samtools stats
    #stats = run_samtools_stats(bam, args.outdir, args.prefix, args.threads)
    print(f"End of workflow for {args.prefix}")

if __name__ == "__main__":
    main()
