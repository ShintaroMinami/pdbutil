from pathlib import Path
from Bio import SeqIO
import dataclasses


@dataclasses.dataclass(frozen=True)
class FastaSequence:
    """
    A dataclass to represent a record in a FASTA file.
    
    Attributes:
        defline (str): The definition line of the sequence.
        sequence (str): The nucleotide or protein sequence.
    """
    defline: list[str]
    sequence: list[str]

    def __iter__(self):
        """
        Returns an iterator over the defline and sequence.
        """
        return iter((self.defline, self.sequence))

    def __len__(self):
        """
        Returns the number of sequences.
        """
        return len(self.sequence)


def read_fasta(file_path: str) -> FastaSequence:
    """
    Reads a .fasta file and returns a list of sequences.

    Args:
        file_path (str): Path to the .fasta file.

    Returns:
        FastaSequence: An object containing lists of deflines and sequences
        extracted from the FASTA file.
    """
    if not Path(file_path).exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    deflines, sequences = [], []
    try:
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                deflines.append(str(record.id))
                sequences.append(str(record.seq))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return FastaSequence(defline=deflines, sequence=sequences)

