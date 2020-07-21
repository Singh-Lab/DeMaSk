#!/usr/bin/env python3

"""
profiles.py
===========

Functions for computing sequence profiles.

"""

import numpy as np

from .utils import AAs


def get_weights(alns: np.ndarray, ident_cutoff: float = 0.8) -> np.ndarray:
    """Get weights for aligned AAs.

    Weights are produced for individual aligned AAs, but currently all
    weights for a given aligned sequence are the same.  The weight for
    a sequence is proportional to the inverse of the number of
    sequences (including itself) above a percent identity cutoff.
    Similar to clustering, this normalizes for uneven phylogenetic
    depths in the database.

    Args:
        alns: An array with dimensions (query length, nseqs + 1)
          containing query + aligned AAs and gaps.
        ident_cutoff: Sequences will be weighted relative to the
          inverse of the number of sequences (including itself) with
          at least this percent identity.

    Returns:
        An array with same dimensions as `alns` containing
        corresponding weights.

    """
    weights = np.zeros(alns.shape)
    residue = np.isin(alns, AAs)
    good = np.empty(alns.shape, dtype=bool)
    AA_match = np.empty(alns.shape, dtype=bool)
    for i in range(alns.shape[1]):
        if alns.shape[1] - 1 > 500 and (i % (alns.shape[1] // 10)) == 0:
            print(".", end="", flush=True)
        good[:] = np.logical_or(residue[:, [i]], residue)
        AA_match[:] = np.logical_and(alns[:, [i]] == alns, good)
        p_idents = np.sum(AA_match, axis=0) / np.sum(good, axis=0)
        p_idents = np.nan_to_num(p_idents)  # Replace nan with 0.
        weights[:, i] = 1.0 / sum(p_idents >= ident_cutoff)
    # Make weights per position (row) sum to 1:
    weights = weights / np.sum(weights, axis=1)[:, np.newaxis]
    return weights


def get_aligned_profile(
    alns: np.ndarray, weights: np.ndarray = None, pseudocount=1e-4
) -> np.ndarray:
    """Create a sequence profile from alignments and optional weights.

    Weights cannot be negative but need not sum to one per position,
    since the profiles will be renormalized per position.  All
    characters other than the 20 canonical amino acids will be
    ignored, so the 20 fractions will always sum to one.

    Args:
        alns: An array with dimensions (query length, num. homologs)
          containing query + aligned AAs and gaps.
        weights: An optional array with same dimensions as `alns`
          containing corresponding weights.
        pseudocount: A small value added to all frequencies before
          normalizing to proportions, allowing the proportions to later
          be log-transformed.

    Returns:
        An array with dimensions (query length, 20) containing the
        proportions of AAs per position.

    """
    if weights is None:
        weights = np.ones(alns.shape)
    else:
        assert np.all(weights >= 0)
    prof = np.zeros((alns.shape[0], len(AAs)))
    for i, aa in enumerate(AAs):
        prof[:, i] = np.sum(np.where(alns == aa, weights, 0), axis=1)
    prof += pseudocount
    prof /= np.sum(prof, axis=1, keepdims=True)
    return prof
