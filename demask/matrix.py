#!/usr/bin/env python3

"""
matrix.py
=========
Code for computing a substitution matrix from DMS data.  The matrix can
be used for variant impact prediction in place of the provided matrix.

"""

import configargparse
import os
import numpy as np

from .utils import AAs, get_filenames, load_dataset


def get_subs_matrix(data: list) -> dict:
    """Compute a substitution matrix from DMS data.

    Values will be averaged by WT + variant identity.  Diagonal
    (synonymous) matrix elements will be set to 0 since the matrix
    represents expected substitution impact.

    Args:
        data: The output of ``load_dataset()``: A list of lists, each containing the residue position,
          WT residue, variant residue, and score for a DMS
          measurement. Position information is not used.

    Returns:
        A dictionary in which keys are (AA1, AA2) tuples.

    """
    # Aggregate scores by WT + variant AA and get averages.
    subs = {(aa1, aa2): [] for aa1 in AAs for aa2 in AAs if aa1 != aa2}
    for d in data:
        subs[d[1], d[2]].append(d[3])
    missing = [sub for sub in subs if len(subs[sub]) == 0]
    assert len(missing) == 0, "Missing substitutions: {}".format(missing)
    matrix = {sub: np.mean(subs[sub]) for sub in subs}
    # Add diagonal elements.
    for aa in AAs:
        matrix[aa, aa] = 0
    return matrix


def write_matrix(matrix: dict, fname: str):
    """Write a substitution matrix to file.

    Args:
        matrix: A dictionary in which keys are (AA1, AA2) tuples, as provided by ``get_subs_matrix()``.
        fname: The file name to write to.

    """
    with open(fname, "w") as out:
        out.write("\t".join(AAs) + "\n")
        for aa1 in AAs:
            vals = ["%.4f" % matrix[aa1, aa2] for aa2 in AAs]
            out.write("\t".join(vals) + "\n")


def prepare_matrix(
    datadir: str,
    outfile: str,
    columns: list = ["pos", "WT", "var", "score"],
):
    """Compute a directional substitution matrix from DMS data.

    Loads a collection of DMS datasets and averages scores by WT +
    variant identity.  Data should already be normalized as desired.
    For the default DeMaSk matrix, variant scores were rank-normalized
    per protein, and then the wild-type fitness levels subtracted, so
    that scores represent delta fitness.

    Dataset files must be tab-delimited with column names.  Columns
    containing residue position, WT residue, variant residue, and
    score must be present in any order, and other columns will be
    ignored.  Entries for which the WT and variant residues are not
    among the 20 canonical amino acids, or are the same, are ignored.

    Args:
        datadir: Name of directory containing (only) DMS data files,
          one per protein.  If a file name is supplied, only that file
          will be used.
        outfile: Name of the file to which the matrix will be written.
        columns: An optional list of column names to read in place of
          'pos', 'WT', 'var', and 'score', respectively.  Must apply
          to all files provided.

    """
    dms = []
    for fname in get_filenames(datadir):
        dms.extend(load_dataset(fname, cols=columns))
    matrix = get_subs_matrix(dms)
    write_matrix(matrix, outfile)


def parse_args():
    p = configargparse.ArgParser(
        description=(
            "Compute a directional amino acid substitution matrix from DMS datasets."
        ),
        epilog="""
        Loads a collection of DMS datasets and averages scores by WT +
        variant identity.  Data should already be normalized as
        desired.  For the default DeMaSk matrix, variant scores were
        rank-normalized per protein, and then the wild-type fitness
        levels subtracted, so that scores represent delta fitness.

        Dataset files must be tab-delimited with column names.
        Columns containing residue position, WT residue, variant
        residue, and score must be present in any order, and other
        columns will be ignored.  Entries for which the WT and variant
        residues are not among the 20 canonical amino acids or are the
        same are ignored.
        """,
    )
    p.add(
        "-d",
        "--datadir",
        required=True,
        help=(
            "Name of directory containing (only) DMS data files, one per "
            "protein.  If a file name is supplied, only that file will be used."
        ),
    )
    p.add("-o", "--outfile", required=True, help="Name of file to write matrix to.")
    p.add(
        "--columns",
        required=False,
        default="pos,WT,var,score",
        help=(
            "An optional comma-separated list of column names to read in place of "
            "'pos', 'WT', 'var', and 'score', respectively.  For example, "
            "'Position,wt,mutant,Fitness'.  Must apply to all files provided."
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
        help=("Configuration file.  Defaults to 'config.ini' in the DeMaSk directory."),
    )
    return p.parse_known_args()[0]


if __name__ == "__main__":
    args = parse_args()
    del args.config
    args.columns = args.columns.split(",")
    prepare_matrix(**vars(args))
