"""Explore oligomer space."""

import logging
import pathlib

import cgexplore as cgx
import matplotlib.pyplot as plt
import numpy as np
import stk
from matplotlib.patches import Patch

import makecages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:  # noqa: C901, D103, PLR0912, PLR0915
    project_dirs = makecages.project.Prj.init_from_json(
        pathlib.Path(__file__).resolve().parent / "project.json"
    )

    conformer_dir = project_dirs.directories["calculation_dir"] / "l1_confs"
    graph_dir = project_dirs.directories["figure_dir"] / "oligo_graphs"
    graph_dir.mkdir(exist_ok=True)

    cis_building_block = stk.BuildingBlock.init_from_file(
        conformer_dir / "c63.mol",
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6H1X3]~[#7X2]~[#6H1X3]",
            bonders=(1,),
            deleters=(),
        ),
    )
    _pd_bb = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=[[0, 0, 0]],
    )
    _capper_bb = stk.BuildingBlock(
        smiles="CC#N",
        functional_groups=(
            stk.SmartsFunctionalGroupFactory(
                smarts="[#7]~[#6]",
                bonders=(0,),
                deleters=(),
            ),
        ),
    )

    graph_details = {}
    for num_pd in (1, 2, 3, 4):
        # Define the iterator.
        for num_cis in range(1, (num_pd * 2) + 1):
            num_capper = ((num_pd * 2) - num_cis) * 2
            if num_capper == 0:
                bb_dict = {
                    _pd_bb: num_pd,
                    cis_building_block: num_cis,
                }

            else:
                bb_dict = {
                    _pd_bb: num_pd,
                    cis_building_block: num_cis,
                    _capper_bb: num_capper,
                }

            iterator = cgx.scram.TopologyIterator(
                building_block_counts=bb_dict
            )
            for topology_code in iterator.yield_graphs():
                # Filter graphs for 1-loops.
                if topology_code.contains_parallels():
                    continue

                name = f"oli_{num_pd}_{num_cis}_{topology_code.idx}"

                graph_details[name] = {
                    "num_cis": num_cis,
                    "t_idx": topology_code.idx,
                    "num_pd": num_pd,
                    "num_monodentate": num_capper,
                    "count_mmmd": 0,
                    "count_mmdd": 0,
                    "count_mddd": 0,
                    "count_dddd": 0,
                }

                makecages.cgexplore.summarise_topology_code(
                    topology_code=topology_code,
                    name=name,
                    figure_dir=graph_dir,
                )

                graph_details[name]["count_doubles"] = (
                    makecages.cgexplore.count_x(topology_code, 4)
                )

                g = topology_code.get_nx_graph()
                n_types = {n: len(nbrs) for n, nbrs in g.adj.items()}

                for nbrs in g.adj.values():
                    types = [n_types[i] for i in list(nbrs.keys())]
                    if len(types) != 4:  # noqa: PLR2004
                        continue

                    count_m = len([i for i in types if i == 1])

                    if count_m == 0:
                        graph_details[name]["count_dddd"] += 1
                    elif count_m == 1:
                        graph_details[name]["count_mddd"] += 1
                    elif count_m == 2:  # noqa: PLR2004
                        graph_details[name]["count_mmdd"] += 1
                    elif count_m == 3:  # noqa: PLR2004
                        graph_details[name]["count_mmmd"] += 1

    logger.info(
        "there are %s graphs in total (ignoring 1-loops!)",
        len(graph_details),
    )

    fig, (ax, ax1, ax2, ax3) = plt.subplots(nrows=4, figsize=(8, 10))

    cmap = {
        1: "tab:gray",
        2: "tab:purple",
        3: "tab:green",
        4: "tab:red",
        5: "tab:cyan",
        6: "tab:pink",
        7: "tab:orange",
        8: "tab:blue",
    }

    alphas = {
        "count_mmmd": 1.0,
        "count_mmdd": 0.7,
        "count_mddd": 0.3,
        "count_dddd": 0.1,
    }

    counts = {}
    counts_pd = {}
    counts_m = {}
    graphs = {}
    counts_alld = 0
    counts_alld_onem = 0
    for cdata in graph_details.values():
        graphs[len(graphs)] = {
            "num_cis": cdata["num_cis"],
            "num_pd": cdata["num_pd"],
            "num_monodentate": cdata["num_monodentate"],
            "count_mddd": cdata["count_mddd"],
            "count_mmmd": cdata["count_mmmd"],
            "count_mmdd": cdata["count_mmdd"],
            "count_dddd": cdata["count_dddd"],
        }

        if cdata["num_cis"] not in counts:
            counts[cdata["num_cis"]] = 0
        counts[cdata["num_cis"]] += 1

        if cdata["num_pd"] not in counts_pd:
            counts_pd[cdata["num_pd"]] = 0
        counts_pd[cdata["num_pd"]] += 1

        if cdata["num_monodentate"] not in counts_m:
            counts_m[cdata["num_monodentate"]] = 0
        counts_m[cdata["num_monodentate"]] += 1

        if (
            cdata["count_mddd"] != 0
            and sum((cdata["count_mmmd"], cdata["count_mmdd"])) == 0
        ):
            counts_alld_onem += 1

        if (
            sum(
                (cdata["count_mddd"], cdata["count_mmmd"], cdata["count_mmdd"])
            )
            == 0
        ):
            counts_alld += 1

    logger.info(
        "count of graphs with all Pd-ditopic %s of %s",
        counts_alld,
        len(graph_details),
    )
    logger.info(
        "count of graphs with all Pd-ditopic and only one Pd-monodentate "
        "%s of %s",
        counts_alld_onem,
        len(graph_details),
    )
    ax2.bar(counts_pd.keys(), counts_pd.values(), ec="k", fc="tab:blue")
    ax3.bar(counts_m.keys(), counts_m.values(), ec="k", fc="tab:blue")

    for nc, c in cmap.items():
        ax.bar(nc, counts[nc], ec="k", fc=c, label=f"L={nc}")

        bottom = np.zeros(
            len([i for i in graphs if graphs[i]["num_cis"] == nc])
        )
        for counter in (
            "count_mmmd",
            "count_mmdd",
            "count_mddd",
            "count_dddd",
        ):
            array = np.array(
                [
                    graphs[i][counter]
                    for i in graphs
                    if graphs[i]["num_cis"] == nc
                ]
            )

            ax1.bar(
                [i for i in graphs if graphs[i]["num_cis"] == nc],
                array,
                ec="k",
                fc=c,
                bottom=bottom,
                alpha=alphas[counter],
            )
            bottom += array

    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel("num. L ($y$)", fontsize=16)
    ax.set_ylabel("num. graphs", fontsize=16)
    ax.legend(fontsize=16, ncol=2)

    ax2.tick_params(axis="both", which="major", labelsize=16)
    ax2.set_xlabel("num. Pd ($x$)", fontsize=16)
    ax2.set_ylabel("num. graphs", fontsize=16)

    ax3.tick_params(axis="both", which="major", labelsize=16)
    ax3.set_xlabel("num. monodentate ($z$)", fontsize=16)
    ax3.set_ylabel("num. graphs", fontsize=16)

    ax1.tick_params(axis="both", which="major", labelsize=16)
    ax1.set_xlabel("graph index", fontsize=16)
    ax1.set_ylabel("count", fontsize=16)
    legend_elements = [
        Patch(facecolor="k", edgecolor="k", alpha=j, label=i.split("_")[1])
        for i, j in alphas.items()
    ]
    ax1.legend(handles=legend_elements, fontsize=16)

    fig.tight_layout()
    fig.savefig(
        project_dirs.directories["figure_dir"] / "graph_properties.png",
        dpi=360,
        bbox_inches="tight",
    )
    fig.savefig(
        project_dirs.directories["figure_dir"] / "graph_properties.pdf",
        dpi=360,
        bbox_inches="tight",
    )
    plt.close()


if __name__ == "__main__":
    main()
