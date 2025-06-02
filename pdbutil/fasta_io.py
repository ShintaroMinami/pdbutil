from pathlib import Path
from Bio import SeqIO
import dataclasses

@dataclasses.dataclass()
class FastaRecord:
    """
    A dataclass to represent a single entry in a FASTA file.
    
    Attributes:
        defline (str): The definition line of the sequence.
        sequence (str): The nucleotide or protein sequence.
    """
    defline: str
    sequence: str


@dataclasses.dataclass()
class FastaRecords:
    entries: list[FastaRecord]

    def __post_init__(self):
        self.sequences = self._sequences()
        self.deflines = self._deflines()

    def __iter__(self):
        """
        Returns an iterator over the entries in the FASTA sequences.
        """
        return iter(self.entries)
    
    def __len__(self):
        """
        Returns the number of entries in the FASTA sequences.
        """
        return len(self.entries)

    def _sequences(self):
        """
        Returns a list of sequences from the FASTA entries.
        """
        return [entry.sequence for entry in self.entries]
    
    def _deflines(self):
        """
        Returns a list of deflines from the FASTA entries.
        """
        return [entry.defline for entry in self.entries]


def read_fasta(file_path: str) -> FastaRecords:
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
    data_list = []
    try:
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                data_list.append(FastaRecord(defline=str(record.id), sequence=str(record.seq)))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return FastaRecords(entries=data_list)

