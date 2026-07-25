import os.path
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "04_data_processing")))
import sascorer

from rdkit import Chem


class TestCase(unittest.TestCase):

  def test_known_small_molecules(self):
    examples = {
      "CCO": 1.980257,
      "c1ccccc1": 1.000000,
    }
    for smi, expected in examples.items():
      mol = Chem.MolFromSmiles(smi)
      self.assertIsNotNone(mol)
      self.assertAlmostEqual(expected, sascorer.calculateScore(mol), 5)


if __name__ == '__main__':
  unittest.main()
