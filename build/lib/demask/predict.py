#!/usr/bin/env python3

"""
predict.py
==========
Functions for generating variant fitness predictions for a query sequence.

"""

import configargparse
import os
import numpy as np

from .profiles import get_weights, get_aligned_profile
from .utils import read_fasta, read_matrix, index, AAs, seq_matrix, entropies


def demask_score(WT: str, var: str, entropy: float, logf_var: float, matrix: dict, coefs: dict) -> float:
    """Calculate the fitness scores for a list of substitutions.

    Args:
        WT: The wild-type amino acid.
        var: The variant amino acid.
        entropy: The Shannon entropy for `WT`'s position in the homolog profile.
        logf_var: The log2-transformed variant AA proportion.
        matrix: A dictionary in which keys are (AA1, AA2) tuples.
        coefs: A dictionary of coefficients with keys "intercept",
          "entropy", "log2f_var", and "matrix".

    Returns:
        The variant's predicted delta fitness.

    """
    score = coefs["intercept"]
    score += coefs["entropy"] * entropy
    score += coefs["log2f_var"] * logf_var
    score += coefs["matrix"] * matrix[WT, var]
    return score


def read_coefficients(fname: str) -> dict:
    """Read DeMaSk model coefficients from file.

    Args:
        fname: File name.

    Returns:
        A dictionary with keys "intercept", "entropy", "log2f_var", and "matrix", and corresponding
        coefficients as values.
    
    """
    with open(fname, "r") as f:
        lines = f.read().splitlines()
        lines = [line.split("\t") for line in lines]
        return {line[0]: float(line[1]) for line in lines}


def compute_scores(wt: str, entropy: list, log_profile: np.ndarray, matrix: dict, coefs: dict) -> list:
    """Compute fitness predictions for all possible substitutions.

    Args:
        wt: The wild-type protein sequence as a string, list, or array
          of characters.
        entropy: A list the length of the query sequence containing
          per-position entropy from the homolog sequence profile.
        log_profile: A sequence profile 2D array as returned by ``get_aligned_profile()``, but
          subsequently log2-transformed.
        matrix: A dictionary in which keys are (AA1, AA2) tuples.
        coefs: A dictionary with keys "intercept", "entropy", "log2f_var", and
          "matrix", and corresponding coefficients as values.

    Returns:
        A list of tuples, each containing the position, WT AA, variant
        AA, score, and the three feature values for each prediction.

    """
    wt = np.array(list(wt))
    posns = index(np.isin(wt, AAs))
    wt = wt[posns]
    entropy = np.array(entropy)[posns]
    log_profile = log_profile[posns, :]

    # Compute scores for all substitutions at all positions.
    variants = [(p, var) for p in range(len(wt)) for var in AAs if wt[p] != var]
    predictions = []
    for i in range(len(variants)):
        p, v = variants[i]
        log2f_var = log_profile[p, index(AAs == v)[0]]
        score = demask_score(wt[p], v, entropy[p], log2f_var, matrix, coefs)
        predictions.append((posns[p] + 1, wt[p], v, score, entropy[p], log2f_var, matrix[wt[p], v]))
    return predictions


def run_demask(
    infile: str,
    matrix: str,
    coefs: str,
    nseqs: int = 500,
    weight_threshold: float = None,
) -> list:
    """Predict fitness for all substitutions in a query sequence.

    Args:
        infile: Name of the file containing a sequence alignment in
          A2M (FASTA) format, with the query protein as the first
          sequence.
        matrix: Name of the file containing the directional
          substitution matrix.
        coefs: Name of the file containing the intercept, entropy, log2f_var, and
          matrix coefficients, with one name and value, separated by
          a tab, per line.
        nseqs: Maximum number of sequences in the alignment to use.
        weight_threshold: Sequence identity threshold used for
          sequence weighting, e.g. 0.8.  Sequences are weighted by the inverse
          of the number of sequences within this percent identity.  If None
          (default), sequence weighting is not used.

    Returns:
        A list of tuples, each containing the position, WT AA, variant
        AA, score, and the three feature values for each prediction.

    """
    seqs = read_fasta(infile, as_dict=False)
    seqs = seqs[:min(len(seqs), nseqs + 1)]  # +1 for query.
    aligned = seq_matrix(seqs)
    query = aligned[:, 0].T     # Query seq minus any gaps.
    aligned = aligned[:, 1:]
    if weight_threshold is None:
        profile = get_aligned_profile(aligned)
    else:
        weights = get_weights(aligned, ident_cutoff=weight_threshold)
        profile = get_aligned_profile(aligned, weights)
    entropy = entropies(profile)
    log_profile = np.log2(profile)
    matrix = read_matrix(matrix)
    coefs = read_coefficients(coefs)
    return compute_scores(query, entropy, log_profile, matrix, coefs)


def parse_args():
    p = configargparse.ArgParser(
        description=(
            "Predict fitness impact of all possible substitutions in "
            "a query protein."
        )
    )
    p.add(
        "-i",
        "--infile",
        required=True,
        help=(
            "Name of the file containing a sequence alignment in A2M (FASTA) "
            "format, with the query protein as the first sequence."
        ),
    )
    p.add("-o", "--outfile", default=None, help="Name of new file to write scores to.")
    dirname = os.path.dirname(__file__)
    config = os.path.join(dirname, "..", "config.ini")
    p.add(
        "-c",
        "--config",
        is_config_file=True,
        metavar="FILE",
        default=config,
        help=("Configuration file.  Defaults to 'config.ini' in the demask directory."),
    )
    matrix = os.path.join(dirname, "..", "data", "matrix.txt")
    p.add(
        "-m",
        "--matrix",
        metavar="FILE",
        default=matrix,
        help=(
            "File containing the directional substitution matrix.  "
            "Defaults to the file at 'DeMaSk/data/matrix.txt'."
        ),
    )
    coefficients = os.path.join(dirname, "..", "data", "coefficients.txt")
    p.add(
        "--coefs",
        metavar="FILE",
        default=coefficients,
        help=(
            "File containing the intercept, entropy, log2f_var, and matrix coefficients.  "
            "Defaults to the file at 'DeMaSk/data/coefficients.txt'."
        ),
    )
    p.add(
        "-n",
        "--nseqs",
        type=int,
        default=500,
        help=(
            "Maximum number of supporting sequences per query sequence.  "
            "Defaults to 500."
        ),
    )
    p.add(
        "-w",
        "--weight_threshold",
        type=float,
        default=None,
        help=(
            "Sequence identity threshold used for "
            "sequence weighting, e.g. 0.8.  Sequences are weighted by the inverse "
            "of the number of sequences within this percent identity.  If None "
            "(default), sequence weighting is not used."
        ),
    )
    args, unknown = p.parse_known_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    outfile = args.outfile
    if outfile is None:
        inbase, inext = os.path.splitext(args.infile)
        outext = ".txt" if inext != ".txt" else "_demask.txt"
        outfile = inbase + outext
    del args.config, args.outfile
    predictions = run_demask(**vars(args))
    with open(outfile, "w") as out:
        out.write("\t".join(
            ["pos", "WT", "var", "score", "entropy", "log2f_var", "matrix"]
        ) + "\n")
        for pred in predictions:
            out.write("{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(*pred))
