import networkx as nx

from dataclasses import dataclass
from pathlib import Path

from chemgraph.io import registry

from . import constants


@dataclass
class ChemGraph:
    name: str | None = None
    """Name of the molecule."""
    graph: nx.Graph = nx.Graph()
    """Graph representation of the molecule."""

    # ============================================================= #

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

    # ============================================================= #

    @classmethod
    def from_file(
        cls,
        path_or_file: Path | str | object,
        name: str | None = None,
        fmt: str | Path | None = None,
        **kwargs,
    ) -> ChemGraph:
        """
        Create a Chemgraph instance from a file.

        Args:
        -----
            path_or_file: Path | str | object
                Path to the file or object to read into a Chemgraph.
            name: str | None
                Default: None.
                Name of the ChemGraph instance.
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
            if not isinstance(path_or_file, (str, Path)):
                raise ValueError(
                    "Format must be specified when path_or_file is not a string or Path."
                )
            fmt = Path(path_or_file).suffix.lstrip(".").lower()

        reader = registry.readers.get(fmt)

        if reader is None:
            raise ValueError(f"No reader registered for format {fmt}")

        data = reader(path_or_file, **kwargs)

        if name is not None:
            data["name"] = name

        return cls(**data)

    # ============================================================= #

    def to_file(
        self, path: str | Path | None = None, fmt: str | Path | None = None, **kwargs
    ) -> None:
        """
        Writes a ChemGraph instance to a file.

        Args:
        -----
            path: Path | str
                Path to the file to write the Chemgraph to.
            fmt: str | None
                Default: None.
                Format of the file.
                If the format is None, extension of the path is used as file format.
                Accepted formats: xyz, mol

        Returns:
        --------
            written: fmt | None
                The specified format is returned if path is not specified.
        """
        if fmt is None:
            if path is None:
                raise ValueError("Format must be specified when path is None.")
            fmt = Path(path).suffix.lstrip(".").lower()

        writer = registry.writers.get(fmt)

        if writer is None:
            raise ValueError(f"No reader registered for format {fmt}")

        if path is not None:
            writer(self, path, **kwargs)
            written = None
        else:
            written = writer(self, **kwargs)

        return written
