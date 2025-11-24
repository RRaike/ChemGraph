import networkx as nx

from dataclasses import dataclass
from pathlib import Path

from chemgraph.io import registry

from . import constants


@dataclass
class ChemGraph:
    name: str | None = None
    """Name of the molecule."""
    graph: nx.Graph | None = nx.Graph()
    """Graph representation of the molecule."""

    def __post_init__(self):
        """
            Normalize graph after initialization:
            - enforce node schema
            - enforce edge schema
            - enforce graph metadata schema
        """
        # === Enforce graph-level schema === #
        for key, default in constants.GRAPH_SCHEMA.items():
            self.graph.graph.setdefault(key, default)

        # Optionally propagate the Chemgraph "name"
        if self.name is not None:
            self.graph.graph.setdefault("name", self.name)

        # === Enforce node schema === #
        for node, attrs in self.graph.nodes(data=True):
            for key, default in constants.NODE_SCHEMA.items():
                attrs.setdefault(key, default)

        # === Enforce edge schema === #
        for u, v, attrs in self.graph.edges(data=True):
            for key, default in constants.EDGE_SCHEMA.items():
                attrs.setdefault(key, default)

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