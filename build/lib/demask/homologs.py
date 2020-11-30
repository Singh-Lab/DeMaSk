#!/usr/bin/env python3

"""
homologs.py
===========
Functions for generating homolog alignments.

"""

import os
import sys
import subprocess
import json
import numpy as np
import configargparse

from .utils import read_fasta, write_fasta


def run_blastp(
    seqfile: str,
    blastp: str,
    db: str,
    evalue: float,
    nseqs: int,
    outfile: str,
    threads: int = 1,
):
    """Execute blastp.

    Executes blastp with the given parameterss and writes output in
    json format (using flag '-outfmt 15').  Throws an exception if the
    blastp run fails.

    Args:
        seqfile: Name of fasta file with query sequence/s.
        blastp: The full path of your blastp executable.
        db: Path of the BLAST database in which supporting sequences
          will be searched.  Must be formatted using the makeblastdb
          program.
        nseqs: Maximum top hits to include per query sequence.
        outfile: Name of results file.
        threads: Number of threads (CPUs) to use.  Defaults to 1.  If 0, it will use a remote database,
          so ``db`` should be the name of an NCBI protein database, e.g. 'nr'.

    """
    command = [
        blastp,
        "-query=" + seqfile,
        "-db=" + db,
        "-num_threads=" + str(threads),
        "-max_target_seqs=" + str(nseqs),
        "-evalue=" + str(evalue),
        "-out=" + outfile,
        "-outfmt=15",
        "-task=blastp-fast", # Comment out when using psiblast
        # "-num_iterations=2" # Only include when using psiblast
    ]
    if threads == 0:
        command[3] = "-remote"
    blrun = subprocess.run(command, env={"BATCH_SIZE": "4000"})
    assert blrun.returncode == 0


def process_query(query: dict, bitscore_cutoff: float = None) -> list:
    """Extract alignment of hits for one query from BLAST output.

    Args:
        query: A dictionary obtained by parsing JSON-formatted BLAST
          output into a nested Python structure and selecting one
          element of the 'BlastOutput2' list.
        bitscore_cutoff: Bitscore (bit_score / query_len) threshold
          for blastp search.  If None (default), no threshold.

    Returns:
        A list of dictionaries, each containing a hit's 'id' and 'seq'.

    """
    hits = query["report"]["results"]["search"]["hits"]
    qlen = query["report"]["results"]["search"]["query_len"]
    # # For psiblast output:
    # hits = query["report"]["results"]["iterations"][-1]["search"]["hits"]
    # qlen = query["report"]["results"]["iterations"][-1]["search"]["query_len"]
    alns = []
    for hit in hits:
        name = hit["description"][0]["id"]
        aln = ["."] * qlen
        pos = hit["hsps"][0]["query_from"] - 1
        qseq = hit["hsps"][0]["qseq"]
        hseq = hit["hsps"][0]["hseq"]
        bitscore = hit["hsps"][0]["bit_score"]
        if bitscore_cutoff is None or bitscore / qlen >= bitscore_cutoff:
            for i in range(len(qseq)):
                if qseq[i] != "-":
                    aln[pos] = hseq[i]
                    pos += 1
            alns.append(dict(id=name, seq="".join(aln)))
    if len(alns) < 10:
        if bitscore_cutoff is not None:
            print(
                "Very few aligned sequences at this bitscore."
                "Trying without bitscore cutoff."
            )
            alns = process_query(query)
        # else:
        #     sys.exit("Error: Not enough aligned sequences.")
    return alns


def get_aligned_AAs(infile: str, bitscore_cutoff: float = None) -> dict:
    """Extract aligned amino acids from BLAST output.

    Args:
        infile: Name of blastp output file in json format (made using
          blastp argument '-outfmt 15').
        bitscore_cutoff: Bitscore (bit_score / query_len) threshold
          for blastp search.  If None (default), no threshold.

    Returns:
        A dictionary in which the keys are the query sequence names
        and the values are lists of aligned hits.  Each hit in the
        list is a dictionary with the hit's 'id' and 'seq'.

    """
    with open(infile, "r") as f:
        blast = json.load(f)["BlastOutput2"]
    aligned = {}
    for query in blast:
        title = query["report"]["results"]["search"]["query_title"]
        # # For psiblast output:
        # title = query["report"]["results"]["iterations"][-1]["search"]["query_title"]
        name = title.split()[0]
        alns = process_query(query, bitscore_cutoff)
        aligned[name] = alns
        if len(alns) == 0:
            print("Warning: No homologs found for {}.".format(name))
        elif len(alns) < 10:
            print("Warning: Very few homologs found for {}.".format(name))
    return aligned


def passes_filters(homolog: str, query: str) -> bool:
    """Check whether a homolog has >= 20% identity with query and whether the alignment is < 10% gaps.

    Args:
        homolog: The aligned homolog sequence. It should be the same length as the query sequence,
          including '.' for end gaps and '-' for internal gaps.
        query: The query sequence.

    Returns:
        True if both conditions are met, otherwise false.

    """
    homolog = np.array(list(homolog))
    query = np.array(list(query))
    identity = np.sum(query == homolog) / np.sum(homolog != ".")
    if identity < 0.2:
        return False
    gap_fraction = np.sum(homolog == "-") / np.sum(homolog != ".")
    if gap_fraction >= 0.1:
        return False
    return True


def find_homologs(
    seqfile: str,
    blastp: str,
    db: str,
    threads: int = 1,
    evalue_cutoff: float = 1e-5,
    bitscore_cutoff: float = None,
    nseqs: int = 500,
    outfile: str = None,
    outdir: str = None,
):
    """Find protein homologs using blastp.

    For each query sequence, an output file in a2m (FASTA) format will
    be produced containing the query sequence and the most similar
    homologs up to a specified number (default 500).

    Args:
        seqfile: Name of fasta file with query sequence/s.
        blastp: The full path of your blastp executable.
        db: Path of the BLAST database in which supporting sequences
          will be searched.  Must be formatted using the makeblastdb
          program.
        threads: Number of threads (CPUs) blastp will use.  Defaults
          to 1.
        evalue_cutoff: E-value threshold for blastp search.  Defaults
          to 1e-5.
        bitscore_cutoff: Bitscore (bit_score / query_len) threshold
          for blastp search.  If None (default), no threshold.
        nseqs: Number of top hits to include per query sequence.
        outfile: Name of output file (only used for single-query
          runs).
        outdir: Name of directory to write output files.  The sequence
          titles (after the ">" and before any whitespace) will be
          used for file names.

    """
    seqbase = os.path.splitext(seqfile)[0]
    blastfile = seqbase + ".blast.json"
    run_blastp(seqfile, blastp, db, evalue_cutoff, nseqs * 2, blastfile, threads)
    alns = get_aligned_AAs(blastfile, bitscore_cutoff)
    # Add query sequence as first in each alignment:
    seqs = read_fasta(seqfile)
    for name in alns:
        # original_count = len(alns[name])
        alns[name] = [x for x in alns[name] if passes_filters(x["seq"], seqs[name])]
        # if len(alns[name]) < original_count:
        #     print("{}: filtered homologs from {} to {}.".format(name, original_count, len(alns[name])))
        alns[name] = alns[name][:min(len(alns[name]), nseqs)]
        alns[name].insert(0, dict(id=name, seq=seqs[name]))

    if outdir is not None:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        for name in alns:
            out = os.path.join(outdir, name + ".a2m")
            write_fasta(alns[name], out)
    else:
        assert len(alns) == 1, "For multiple queries, provide outdir."
        if outfile is None:
            outfile = seqbase + ".a2m"
        write_fasta(list(alns.values())[0], outfile)


def parse_args():
    p = configargparse.ArgParser(
        description="Find protein homologs using blastp.",
        epilog="""
        For each query sequence, an output file in a2m (FASTA) format will
        be produced containing the query sequence and the most similar
        homologs up to a specified number (default 500).
        """,
    )
    p.add(
        "-s",
        "--seqfile",
        required=True,
        help=(
            "Name of a fasta file with one or more query sequences."
        ),
    )
    dirname = os.path.dirname(__file__)
    config = os.path.join(dirname, "..", "config.ini")
    p.add(
        "-c",
        "--config",
        is_config_file=True,
        metavar="FILE",
        default=config,
        help=(
            "Configuration file (e.g. the same configuration file used for "
            "running DeMaSk).  Defaults to 'config.ini' in the demask directory."
        ),
    )
    p.add("--blastp", required=True, help=("The full path of your blastp executable."))
    p.add(
        "--db",
        required=True,
        metavar="PATH",
        help=(
            "Path of the BLAST database in which supporting sequences will "
            "be searched.  Must be formatted using the makeblastdb program."
        ),
    )
    p.add(
        "-t",
        "--threads",
        type=int,
        default=1,
        help=(
            "Number of threads (CPUs) blastp will use.  Defaults to 1."
            "If 0, it will use a remote database, "
            "so 'db' should be the name of an NCBI protein database, e.g. 'nr'."
        ),
    )
    p.add(
        "-e",
        "--evalue_cutoff",
        type=float,
        default=1e-5,
        help="E-value threshold for blastp hits.  Defaults to 1e-5.",
    )
    p.add(
        "-b",
        "--bitscore_cutoff",
        type=float,
        default=None,
        help="Bits per query residue threshold for blastp hits.  Default is no threshold",
    )
    p.add(
        "-n",
        "--nseqs",
        type=int,
        default=500,
        help=(
            "Maximum number of top hits to include per "
            "query sequence.  Defaults to 500."
        ),
    )
    out = p.add_mutually_exclusive_group(required=True)
    out.add(
        "-o",
        "--outfile",
        help=("Name of output file (only used for single-query runs)."),
    )
    out.add(
        "-d",
        "--outdir",
        help=(
            "Name of directory to write output files in a2m (FASTA) format.  "
            "The sequence titles (after the '>' and before any whitespace) "
            "will be used for file names."
        ),
    )
    return p.parse_known_args()[0]


if __name__ == "__main__":
    args = parse_args()
    del args.config
    find_homologs(**vars(args))
