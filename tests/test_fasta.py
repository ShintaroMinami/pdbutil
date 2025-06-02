import unittest
from pathlib import Path
try:
    from pdbutil.fasta_io import FastaSequence, read_fasta
except ImportError:
    import sys
    sys.path.append(str(Path(__file__).resolve().parent.parent))
    from pdbutil.fasta_io import FastaSequence, read_fasta
dir_path = Path(__file__).resolve().parent

test_files = {
    "test_fasta": dir_path / "data/test.fasta"
    }

class TestFastaSequence(unittest.TestCase):

    def test_iter(self):
        """Test the __iter__ method."""
        fasta = FastaSequence(defline=["seq1"], sequence=["SEQUENCE"])
        self.assertEqual(list(iter(fasta)), [["seq1"], ["SEQUENCE"]])

    def test_len(self):
        """Test the __len__ method."""
        fasta = FastaSequence(defline=["seq1", "seq2"], sequence=["SEQUENCEONE", "SEQUENCETWO"])
        self.assertEqual(len(fasta), 2)


class TestReadFasta(unittest.TestCase):

    def test_read_fasta(self):
        """Test the read_fasta function with a valid FASTA file."""
        test_file = test_files["test_fasta"]
        with open(test_file, "w") as f:
            f.write(">seq1\nSEQUENCEONE\n>seq2\nSEQUENCETWO\n")

        fasta = read_fasta(test_file)

        self.assertEqual(len(fasta), 2)
        self.assertEqual(fasta.defline, ["seq1", "seq2"])
        self.assertEqual(fasta.sequence, ["SEQUENCEONE", "SEQUENCETWO"])

        Path(test_file).unlink()

    def test_read_fasta_file_not_found(self):
        """Test the read_fasta function with a non-existent file."""
        with self.assertRaises(FileNotFoundError):
            read_fasta("non_existent.fasta")

if __name__ == "__main__":
    unittest.main()
