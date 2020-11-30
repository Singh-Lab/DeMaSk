"""
utils.py
========
Utility functions for the demask package.

"""

import os
import numpy as np


AAs = np.array(list("ACDEFGHIKLMNPQRSTVWY"))


def read_fasta(filename: str, as_dict: bool = True) -> dict:
    """Read one or more sequences from a fasta file.

    Names are extracted from the header line after '>' and before any
    whitespace.

    Args:
        filename: Name of the fasta file to read.
        as_dict: If set to false, will return a list of sequences.

    Returns:
        A dictionary with names as keys and sequence strings as
        values, or a list of sequences.

    """
    seqs = []
    with open(filename, "r") as f:
        lines = f.read().splitlines()
    lines = [line for line in lines if line != ""]
    name = None
    for line in lines:
        if line[0] == ">":
            name = line[1:].split()[0]
            seqs.append(dict(name=name, seq=""))
        else:
            seqs[-1]["seq"] += line.upper()
    if as_dict:
        seqs = {x["name"]: x["seq"] for x in seqs}
    else:
        seqs = [x["seq"] for x in seqs]
    return seqs


def write_fasta(seqs: list, filename: str):
    """Write one or more sequences to a fasta file.

    Args:
        seqs: A list of dictionaries, each containing a sequence's
          'id' and 'seq'.
        filename: Name of fasta file to write to.

    """
    with open(filename, "w") as out:
        for seq in seqs:
            out.write(">" + seq["id"] + "\n" + seq["seq"] + "\n")


def seq_matrix(seqs: list) -> np.ndarray:
    """Get a residue matrix from a list of sequences.
    
    Returns a numpy array of characters with dimensions (seq length,
    nseqs).  The rows correspond to letters in the first (query)
    sequence, so positions in which the first sequence has a gap ('-'
    or '.') are removed.  Lowercase (masked) residues are capitalized
    so that they are counted as valid AAs.

    """
    alns = np.full((len(seqs[0]), len(seqs)), ".", dtype="U1")
    for i, seq in enumerate(seqs):
        alns[:, i] = list(seq.upper())  # Convert lowercase (masked) AAs to upper.
    alns = alns[np.isin(alns[:, 0], ["-", "."], invert=True)]
    return alns


def extract_DMS_WT(datadir: str, output: str):
    """Infer query sequence/s based on substitutions listed in DMS data.

    Amino acids specified by the 'pos' and 'WT' columns will be
    assembled, and any missing positions will be represented as 'X' in
    the sequence.  Sequences will be written to a fasta file.

    This method is not used by DeMaSk, but is provided for convenience if the WT sequence for
    a DMS dataset is missing.

    Args:
        datadir: Name of the directory containing (only) DMS data
          files, one per protein.  If a file name is supplied, only
          that file will be used.
        output: Name of the output (fasta) file.

    """
    seqs = dict()
    if os.path.isdir(datadir):
        fnames = [os.path.join(datadir, f) for f in os.listdir(datadir)]
    elif os.path.isfile(datadir):
        fnames = [datadir]
    for fname in fnames:
        lines = open(fname, "r").read().splitlines()
        lines = [line.split("\t") for line in lines]
        pos_col = lines[0].index("pos")
        WT_col = lines[0].index("WT")
        pos = [int(line[pos_col]) for line in lines[1:]]
        WT = [line[WT_col] for line in lines[1:]]
        chars = ["X"] * max(pos)
        for i in range(len(pos)):
            chars[pos[i] - 1] = WT[i]
        name = os.path.split(fname)[1]
        name = os.path.splitext(name)[0]
        name = name.replace(" ", "_")  # No spaces in fasta seq IDs.
        seqs[name] = "".join(chars)
    with open(output, "w") as out:
        for name in sorted(seqs.keys()):
            out.write(">" + name + "\n" + seqs[name] + "\n\n")


def read_matrix(fname: str) -> dict:
    """Read a substitution matrix from file.

    The file must be tab-delimited with column names in the first
    line.  Row names should not be provided and will be copied from
    the column names.  If the matrix is not symmetric, rows should
    correspond to AA1 and columns to AA2.

    Args:
        fname: The file name.

    Returns:
        A dictionary in which keys are (AA1, AA2) tuples.

    """
    matrix = dict()
    with open(fname, "r") as f:
        lines = f.read().splitlines()
    mat_AAs = lines[0].split("\t")
    for i, aa1 in enumerate(mat_AAs):
        vals = lines[i + 1].split("\t")
        for j, aa2 in enumerate(mat_AAs):
            matrix[aa1, aa2] = float(vals[j])
    return matrix


def get_filenames(location: str) -> list:
    """Get file names from user-specified directory.

    If the input is actually the name of an existing file, only that
    file name will be returned.

    Args:
        location: The name of the directory or file.

    Returns:
        A list of file names.

    """
    if os.path.isdir(location):
        fnames = [os.path.join(location, f) for f in os.listdir(location) if not f.startswith('.')]
    elif os.path.isfile(location):
        fnames = [location]
    else:
        raise Exception(f"Data file/s not found: {location}")
    return fnames


def load_dataset(fname: str, cols: list = ["pos", "WT", "var", "score"]) -> list:
    """Load a DMS dataset from a file.

    File must be tab-delimited with column names.  Columns containing
    residue position, WT residue, variant residue, and score must be
    present in any order, and other columns will be ignored.  Entries
    for which the WT and variant residues are not among the 20
    canonical amino acids or are the same are removed.

    Args:
        fname: Name of the file.
        cols: An optional list of column names to read in place of
          'pos', 'WT', 'var', and 'score', respectively.

    Returns:
        A list of lists, each containing the residue position, WT
        residue, variant residue, and score for a data point.

    """
    with open(fname, "r") as f:
        lines = f.read().splitlines()
    lines = [line.split("\t") for line in lines]
    pos_col = lines[0].index(cols[0])
    WT_col = lines[0].index(cols[1])
    var_col = lines[0].index(cols[2])
    score_col = lines[0].index(cols[3])
    lines = [line for line in lines[1:] if line[score_col] != "NA"]
    dms = [
        [int(line[pos_col]), line[WT_col], line[var_col], float(line[score_col])]
        for line in lines
    ]
    dms = [d for d in dms if d[1] in AAs and d[2] in AAs and d[1] != d[2]]
    return dms


def index(x: np.ndarray) -> np.ndarray:
    """Get array of indices for True elements in 1D array."""
    assert x.ndim == 1
    return np.where(x)[0]


def entropies(profile: np.ndarray) -> list:
    """Calculate Shannon entropy (log base 2) for each position of a sequence profile.

    Args:
        profile: A 2D array of AA proportions as returned by ``get_aligned_profile()``.

    Returns:
        A list of entropy values corresponding to the rows of ``profile``.

    """
    return [-sum(p[np.nonzero(p)] * np.log2(p[np.nonzero(p)])) for p in profile]
