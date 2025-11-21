import networkx as nx

from dataclasses import dataclass
from pathlib import Path

from chemgraph.io import registry

@dataclass
class Chemgraph:
    name: str | None = None
    """Name of the molecule."""
    graph: nx.Graph | None = None
    """Graph representation of the molecule."""

    @classmethod
    def from_file(cls, path: str | Path,  fmt: str | Path | None = None) -> Chemgraph:
        """
            Create a Chemgraph instance from a .xyz file.

            Args:
            -----
                path: Path | Str
                    Path to the file to read into a Chemgraph.
                fmt: str | None
                    Default: None.
                    Format of the file. 
                    If the format is None, extension of the path is used as file format.
                    Accepted formats: xyz

            Returns:
            --------
                Chemgraph: 
                    New ChemGraph instance with data from the .xyz file. 
        """
        if fmt is None:
            fmt = Path(path).suffix.lstrip('.').lower()

        reader = registry.readers.get(fmt)

        if reader is None:
            raise ValueError("No reader registered for format .xyz")
        data = reader(path)
        return cls(**data)