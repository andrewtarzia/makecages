"""Utilities for cgx topologies."""

import pathlib
from collections import Counter

import cgexplore as cgx
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import rustworkx as rx


def count_x(topology_code: cgx.scram.TopologyCode, x: int) -> dict:
    """Count the number of x loops."""
    weighted_graph = topology_code.get_weighted_graph()

    filtered_paths = set()
    for node in weighted_graph.nodes():
        paths = list(
            rx.graph_all_simple_paths(
                weighted_graph,
                origin=node,  # type: ignore[call-arg]
                to=node,  # type: ignore[call-arg]
                cutoff=12,
                min_depth=4,
            )
        )

        for path in paths:
            # This does not currently handle repeated cycles, inverse
            # cycles and shifted cycles.
            if (
                tuple(path) not in filtered_paths
                and tuple(path[::-1]) not in filtered_paths
            ):
                filtered_paths.add(tuple(path))
    path_lengths = [len(i) - 1 for i in filtered_paths]

    return Counter(path_lengths)[x]


def summarise_topology_code(
    topology_code: cgx.scram.TopologyCode,
    name: str,
    figure_dir: pathlib.Path,
) -> None:
    """Use networkx to layout and summarise a graph."""
    g = topology_code.get_nx_graph()
    degree_sequence = sorted((d for n, d in g.degree()), reverse=True)

    fig = plt.figure(figsize=(8, 5))
    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(1, 2)

    ax0 = fig.add_subplot(axgrid[:, :1])
    gcc = g.subgraph(
        sorted(nx.connected_components(g), key=len, reverse=True)[0]
    )
    pos = nx.spring_layout(gcc, seed=10396953)
    nx.draw_networkx_nodes(gcc, pos, ax=ax0, node_size=20)
    nx.draw_networkx_edges(gcc, pos, ax=ax0, alpha=0.4)
    ax0.tick_params(axis="both", which="major", labelsize=16)
    ax0.set_title(f"Connected components of {name}", fontsize=16)
    ax0.set_axis_off()

    ax2 = fig.add_subplot(axgrid[:, 1:])
    ax2.bar(*np.unique(degree_sequence, return_counts=True))
    ax2.tick_params(axis="both", which="major", labelsize=16)
    ax2.set_title("Degree histogram", fontsize=16)
    ax2.set_xlabel("Degree", fontsize=16)
    ax2.set_ylabel("# of Nodes", fontsize=16)

    fig.tight_layout()
    fig.savefig(
        figure_dir / f"g_{name}.png",
        dpi=360,
        bbox_inches="tight",
    )
    plt.close()
