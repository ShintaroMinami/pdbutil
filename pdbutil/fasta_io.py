from pathlib import Path
from Bio import SeqIO
import dataclasses
from typing import Iterable

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
    """
    A dataclass to represent multiple entries in a FASTA file.

    Attributes:
        entries (list[FastaRecord]): A list of FastaRecord objects.
    """
    entries: list[FastaRecord]

    def __post_init__(self):
        """
        Initializes sequences and deflines attributes after the dataclass is created.
        """
        self.sequences: list[str] = self._sequences()
        self.deflines: list[str] = self._deflines()

    def __iter__(self) -> Iterable[FastaRecord]:
        """
        Returns an iterator over the entries in the FASTA sequences.

        Returns:
            Iterable[FastaRecord]: An iterator over FastaRecord objects.
        """
        return iter(self.entries)

    def __len__(self) -> int:
        """
        Returns the number of entries in the FASTA sequences.

        Returns:
            int: The number of entries.
        """
        return len(self.entries)

    def __str__(self) -> str:
        """
        Returns a string representation of the FASTA records.

        Returns:
            str: A string representation of the FASTA records in FASTA format.
        """
        return "\n".join(f">{entry.defline}\n{entry.sequence}" for entry in self.entries)

    def _sequences(self) -> list[str]:
        """
        Returns a list of sequences from the FASTA entries.

        Returns:
            list[str]: A list of sequences.
        """
        return [entry.sequence for entry in self.entries]
    
    def _deflines(self) -> list[str]:
        """
        Returns a list of deflines from the FASTA entries.

        Returns:
            list[str]: A list of deflines.
        """
        return [entry.defline for entry in self.entries]


def read_fasta(file_path: str) -> FastaRecords:
    """
    Reads a .fasta file and returns a FastaRecords object containing the sequences and deflines.

    Args:
        file_path (str): Path to the .fasta file.

    Returns:
        FastaRecords: An object containing lists of deflines and sequences extracted from the FASTA file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    
    if not Path(file_path).exists():
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    data_list: list[FastaRecord] = []
    
    try:
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                data_list.append(FastaRecord(defline=str(record.id), sequence=str(record.seq)))
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    
    return FastaRecords(entries=data_list)

