#!/usr/bin/env python3

"""
fit.py
======
Functions for fitting a linear model to relate entropy, variant frequency, and a substitution matrix to
DMS scores.

"""

import os
import numpy as np
import configargparse

from .utils import read_fasta, read_matrix, get_filenames, load_dataset, seq_matrix, entropies, index, AAs
from .profiles import get_weights, get_aligned_profile


def coefficients(dms: dict, entropy: dict, log_profile: dict, matrix: dict) -> dict:
    """Compute coefficients for the DeMaSk linear model.

    Args:
        dms: A dictionary in which keys are dataset names and values are DMS data as returned
          by ``load_dataset()``.
        entropy: A dictionary in which keys are dataset names and values are iterables with
          entropy for each protein position.
        log_profile: A dictionary in which keys are dataset names and values are
          the sequence profile 2D arrays returned by ``get_aligned_profile()``.
        matrix: A dictionary in which keys are (AA1, AA2) tuples.

    Returns:
        A dictionary with keys "intercept", "entropy", "log2f_var", and "matrix", and corresponding
        coefficients as values.

    """
    features = []
    scores = []
    for dataset in dms:
        for d in dms[dataset]:
            log2f_var = log_profile[dataset][d[0] - 1, index(AAs == d[2])[0]]
            # Include 1 for intercept term.
            features.append([1, entropy[dataset][d[0] - 1], log2f_var, matrix[d[1], d[2]]])
            scores.append(d[3])
    X = np.matrix(features)
    y = np.matrix(scores).T
    beta = np.linalg.inv(X.T * X) * X.T * y
    beta = np.array(beta.T)[0]
    return dict(intercept=beta[0], entropy=beta[1], log2f_var=beta[2], matrix=beta[3])


def fit_model(
    datadir: str,
    alndir: str,
    matrix: str,
    outfile: str,
    columns: list = ["pos", "WT", "var", "score"],
    nseqs: int = 500,
    weight_threshold: float = None,
):
    """Fit the DeMaSk model coefficients.

    For a set of DMS datasets, load the DMS scores, aligned homologs,
    and pre-computed substitution matrix, and fit a linear model,
    saving the coefficients to a file for later prediction on new
    proteins.

    Dataset files must be tab-delimited with column names.  Columns
    containing residue position, WT residue, variant residue, and
    score must be present in any order, and other columns will be
    ignored.  Entries for which the WT and variant residues are not
    among the 20 canonical amino acids or are the same are ignored.

    DMS scores should already be normalized as desired.  For the
    default DeMaSk matrix, variant scores were rank-normalized per
    protein, and then the wild-type fitness levels subtracted, so that
    scores represent delta fitness.

    Args:
        datadir: Name of the directory containing (only) DMS data
          files, one per protein.
        alndir: Name of the directory containing alignment files, each
          of which must be named as the basename of a DMS data file
          followed by the .a2m extension.
        matrix: Name of the file containing the directional
          substitution matrix.
        outfile: Name of the file in which to save the coefficients.
        columns: An optional list of DMS data column names to read in
          place of 'pos', 'WT', 'var', and 'score', respectively.
          Must apply to all files provided.
        nseqs: Maximum number of sequences in the alignment to use.
        weight_threshold: Sequence identity threshold used for
          sequence weighting, e.g. 0.8.  Sequences are weighted by the inverse
          of the number of sequences within this percent identity.  If None
          (default), sequence weighting is not used.

    """
    dms = {}
    entropy = {}
    log_profile = {}
    for fname in get_filenames(datadir):
        name = os.path.splitext(os.path.basename(fname))[0]
        dms[name] = load_dataset(fname, cols=columns)
        alnfile = os.path.join(alndir, name + ".a2m")
        seqs = read_fasta(alnfile, as_dict=False)
        if nseqs is not None:
            seqs = seqs[: min(len(seqs), nseqs + 1)]  # +1 for query.
        aligned = seq_matrix(seqs)
        aligned = aligned[:, 1:]  # Remove query seq from homologs.
        if weight_threshold is None:
            profile = get_aligned_profile(aligned)
        else:
            weights = get_weights(aligned, ident_cutoff=weight_threshold)
            profile = get_aligned_profile(aligned, weights)
        entropy[name] = entropies(profile)
        log_profile[name] = np.log2(profile)
    matrix = read_matrix(matrix)
    coefs = coefficients(dms, entropy, log_profile, matrix)
    with open(outfile, "w") as f:
        for term in ["intercept", "entropy", "log2f_var", "matrix"]:
            f.write("{}\t{:.6f}\n".format(term, coefs[term]))


def parse_args():
    p = configargparse.ArgParser(
        description=("Compute coefficients for the DeMaSk predictor."),
        epilog="""
        For a set of DMS datasets, load the DMS scores, aligned homologs,
        and pre-computed substitution matrix, and fit a linear model,
        saving the coefficients to a file for later prediction on new
        proteins.
        
        Dataset files must be tab-delimited with column names.
        Columns containing residue position, WT residue, variant
        residue, and score must be present in any order, and other
        columns will be ignored.  Entries for which the WT and variant
        residues are not among the 20 canonical amino acids or are the
        same are ignored.

        DMS scores should already be normalized as desired.  For the
        default DeMaSk matrix, variant scores were rank-normalized per
        protein, and then the wild-type fitness levels subtracted, so that
        scores represent delta fitness.
        """,
    )
    p.add(
        "-d",
        "--datadir",
        required=True,
        help="Name of directory containing (only) DMS data files, one per protein.",
    )
    p.add(
        "-a",
        "--alndir",
        required=True,
        help=(
            "Name of the directory containing alignment files, each"
            "of which must be named as the basename of a DMS data file"
            "followed by the '.a2m' extension."
        ),
    )
    dirname = os.path.dirname(__file__)
    matrix = os.path.join(dirname, "..", "data", "matrix.txt")
    p.add(
        "-m",
        "--matrix",
        metavar="FILE",
        default=matrix,
        help=(
            "File containing the directional substitution matrix.  "
            "Defaults to the file at 'demask/matrix/matrix.txt'."
        ),
    )
    p.add(
        "-o",
        "--outfile",
        required=True,
        help="Name of the file in which to save the coefficients.",
    )
    p.add(
        "--columns",
        required=False,
        default="pos,WT,var,score",
        help=(
            "An optional comma-separated list of DMS data column names to read in place of "
            "'pos', 'WT', 'var', and 'score', respectively.  For example, "
            "'Position,wt,mutant,Fitness'.  Must apply to all files provided."
            ),
    )
    config = os.path.join(dirname, "..", "config.ini")
    p.add(
        "-c",
        "--config",
        is_config_file=True,
        metavar="FILE",
        default=config,
        help="Configuration file.  Defaults to 'config.ini' in the demask directory.",
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
    del args.config
    args.columns = args.columns.split(",")
    fit_model(**vars(args))
