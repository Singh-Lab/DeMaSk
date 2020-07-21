import unittest
import json
import os
import numpy as np

from demask.utils import entropies, load_dataset, read_fasta, read_matrix, seq_matrix
from demask.matrix import get_subs_matrix, prepare_matrix
from demask.homologs import process_query, find_homologs
from demask.profiles import get_aligned_profile, get_weights
from demask.predict import compute_scores, read_coefficients


class TestMatrix(unittest.TestCase):
    def setUp(self):
        if os.path.exists("test/tmp.txt"):
            os.remove("test/tmp.txt")

    def test_load_dataset_col_order(self):
        # Ensure data cols can be in any order.
        d1 = load_dataset("test/data.txt")
        d2 = load_dataset("test/data_order.txt")
        self.assertEqual(d1, d2)

    def test_load_dataset_col_names(self):
        # Ensure custom column names can be specified.
        d1 = load_dataset("test/data.txt")
        d2 = load_dataset(
            "test/data_colnames.txt", cols=["Position", "wt", "Mutant", "fitness"]
        )
        self.assertEqual(d1, d2)

    def test_load_dataset_filtering(self):
        # Ensure synonymous variants and those with non-canonical WT
        # or variant AAs are ignored.
        d = load_dataset("test/data.txt")
        self.assertEqual(len(d), 380)

    def test_get_subs_matrix(self):
        d = load_dataset("test/data.txt")
        mat = get_subs_matrix(d)
        self.assertEqual(len(mat), 400)

    def test_get_subs_matrix_incomplete(self):
        # Ensure understandable exception is raised when a WT/var
        # combo is absent from the data.
        d = load_dataset("test/data_incomplete.txt")
        self.assertRaises(AssertionError, get_subs_matrix, d)

    def test_prepare_matrix_file(self):
        # Ensure matrix file is produced from a single data file.
        prepare_matrix("test/data.txt", "test/tmp.txt")
        self.assertTrue(os.path.isfile("test/tmp.txt"))

    def test_prepare_matrix_dir(self):
        # Ensure matrix file is produced from a directory of data.
        prepare_matrix("test/data", "test/tmp.txt")
        self.assertTrue(os.path.isfile("test/tmp.txt"))

    def tearDown(self):
        if os.path.exists("test/tmp.txt"):
            os.remove("test/tmp.txt")


class TestHomologs(unittest.TestCase):
    def test_process_query(self):
        with open("test/P46937.blast.json", "r") as f:
            query = json.load(f)["BlastOutput2"][0]
        hits = process_query(query)
        self.assertTrue(len(hits) == 50)

    def test_no_blast_hits(self):
        ## Uncomment this to actually run blastp:
        # find_homologs(
        #     "test/random.fa",
        #     blastp="/usr/local/ncbi/blast/bin/blastp",
        #     db="/Volumes/storage/db/uniref90.fasta",
        #     threads=2,
        #     outfile="test/random.a2m",
        # )
        seqs = read_fasta("test/random.a2m")
        self.assertEqual(len(seqs), 1)


class TestProfiles(unittest.TestCase):
    def test_get_weights(self):
        seqs = read_fasta("test/P46937.a2m", as_dict=False)
        aln = seq_matrix(seqs)[:, 1:]
        wts = get_weights(aln)
        self.assertTrue(np.all(wts > 0) and np.all(wts < 1))

    def test_get_aligned_profile(self):
        # Ensure sum of each row of output array is > 0 and <= 1.
        seqs = read_fasta("test/P46937.a2m", as_dict=False)
        aln = seq_matrix(seqs)[:, 1:]
        wts = get_weights(aln)
        profile = get_aligned_profile(aln, wts)
        sums = np.sum(profile, axis=1)
        self.assertTrue(np.all(sums > 0) and np.all(sums <= 1.01))

    def test_gaps_in_query(self):
        # demask should accept gaps in query, e.g. for true MSA input.
        seqs1 = read_fasta("test/test_aln_ungapped.a2m", as_dict=False)
        aln1 = seq_matrix(seqs1)
        seqs2 = read_fasta("test/test_aln_gapped.a2m", as_dict=False)
        aln2 = seq_matrix(seqs2)
        self.assertTrue(aln1.shape == aln2.shape and np.all(aln1 == aln2))

    def test_no_homologs(self):
        seqs = read_fasta("test/random.a2m", as_dict=False)
        aln = seq_matrix(seqs)[:, 1:]
        wts = get_weights(aln)
        profile = get_aligned_profile(aln, wts)
        self.assertTrue(np.all(profile == 1/20))


class TestPredict(unittest.TestCase):
    def setUp(self):
        self.seqs = read_fasta("test/P46937.a2m", as_dict=False)
        aln = seq_matrix(self.seqs)[:, 1:]
        profile = get_aligned_profile(aln)
        self.entropies = entropies(profile)
        self.log_profile = np.log2(profile + 1 / profile.shape[1])
        self.matrix = read_matrix("data/matrix.txt")
        self.coefs = read_coefficients("data/coefficients.txt")

    def test_compute_scores_string(self):
        # Ensure output file is written given WT sequence as a string.
        wt = self.seqs[0]
        predictions = compute_scores(wt, self.entropies, self.log_profile, self.matrix, self.coefs)
        self.assertEqual(len(predictions), len(self.seqs[0]) * 19)

    def test_compute_scores_list(self):
        # Ensure output file is written given WT sequence as a list.
        wt = list(self.seqs[0])
        predictions = compute_scores(wt, self.entropies, self.log_profile, self.matrix, self.coefs)
        self.assertEqual(len(predictions), len(self.seqs[0]) * 19)

    def test_compute_scores_array(self):
        # Ensure output file is written given WT sequence as an array.
        wt = np.array(list(self.seqs[0]))
        predictions = compute_scores(wt, self.entropies, self.log_profile, self.matrix, self.coefs)
        self.assertEqual(len(predictions), len(self.seqs[0]) * 19)


if __name__ == "__main__":
    unittest.main()
