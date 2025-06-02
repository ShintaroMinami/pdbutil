import unittest
from pathlib import Path
import numpy as np
try:
    from pdbutil.rmsd import kabsch, superpose, calc_rmsd
except ImportError:
    import sys
    sys.path.append(str(Path(__file__).resolve().parent.parent))
    from pdbutil.rmsd import kabsch, superpose, calc_rmsd

class TestRMSD(unittest.TestCase):

    def test_kabsch(self):
        """Test the Kabsch algorithm."""
        A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 1]])
        B = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 1, 1]])
        self.assertGreater(len(A), 3, "A must have at least 3 points")
        self.assertGreater(len(B), 3, "B must have at least 3 points")
        X, Y = kabsch(A, B)
        self.assertEqual(X.shape, A.shape)
        self.assertEqual(Y.shape, B.shape)
        rmsd = np.sqrt(((X-Y)**2).sum(axis=-1).mean())
        self.assertAlmostEqual(rmsd, 0.79473, places=5, msg="RMSD should be close to 0.794731876... for this test case")

    def test_superpose(self):
        """Test the superpose function."""
        ref = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 1]])
        trg = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 1, 1]])
        self.assertGreater(len(ref), 3, "Reference must have at least 3 points")
        self.assertGreater(len(trg), 3, "Target must have at least 3 points")
        superposed, rmsd = superpose(ref, trg, return_rmsd=True)
        self.assertEqual(superposed.shape, trg.shape)
        self.assertTrue(isinstance(rmsd, np.ndarray), "RMSD should be a float")

    def test_calc_rmsd(self):
        """Test the calc_rmsd function."""
        A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 1]])
        B = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0], [1, 1, 1]])
        self.assertGreaterEqual(len(A), 3, "A must have at least 3 points")
        self.assertGreaterEqual(len(B), 3, "B must have at least 3 points")
        rmsd = calc_rmsd(A, B)
        self.assertTrue(isinstance(rmsd, np.ndarray))
        self.assertAlmostEqual(rmsd, 0.79473, places=5, msg="RMSD should be close to 0.794731876... for this test case")

        
if __name__ == "__main__":
    unittest.main()
