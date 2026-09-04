"""Build cages."""

import logging
import pathlib
from collections import defaultdict

import atomlite
import cgexplore as cgx
import matplotlib.pyplot as plt
import stk
import stko

import makecages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)

topology_codes = {
    "m2l4a": cgx.scram.get_stk_topology_code(stk.cage.M2L4Lantern),
    "m3l6a": cgx.scram.get_stk_topology_code(stk.cage.M3L6),
}


def plotting(project_dirs: makecages.project.Prj) -> None:  # noqa: C901, PLR0912, PLR0915
    """Make some plots."""
    posmaps = {
        "m2l4": {"cis": 0, "trans": 1},
        "m3l6": {
            "cis": 0,
            "trans": 1,
            "trans2": 2,
            "trans3": 3,
        },
    }
    mixposmaps = {
        "m2l4": {},
        "m3l6": {
            "trans2": 1,
            "trans2-4": 2,
            "trans2-6": 3,
            "trans2-9": 4,
            "trans2-10": 5,
            "trans2-11": 6,
            "trans": 7,
            "trans2-5": 0,
            "trans2-12": -1,
            "trans2-13": -2,
            "trans2-14": -3,
            "trans2-15": -4,
            "trans3": -5,
        },
    }
    vaposmaps = {
        "m2l4": {},
        "m3l6": {
            "trans2": 0,
            "trans2-1": 1,
            "trans2-1a": 2,
            "trans2-2": 3,
            "trans2-3": 4,
            "trans2-7": 5,
            "trans2-8": 6,
        },
    }
    name_map = {
        "trans": "$t$-1",
        "trans2": "$t$-2",
        "trans3": "$t$-3",
        "cis": "$c$",
        "trans2-1": "$t$-2-1",
        "trans2-1a": "$t$-2-1a",
        "trans2-2": "$t$-2-2",
        "trans2-3": "$t$-2-3",
        "trans2-7": "$t$-2-7",
        "trans2-8": "$t$-2-8",
        "trans2-4": "5($t$-2)+1($t$-1)",
        "trans2-6": "4($t$-2)+2($t$-1)",
        "trans2-9": "3($t$-2)+3($t$-1)",
        "trans2-10": "2($t$-2)+4($t$-1)",
        "trans2-11": "1($t$-2)+5($t$-1)",
        "trans2-5": "5($t$-2)+1($t$-3)",
        "trans2-12": "4($t$-2)+2($t$-3)",
        "trans2-13": "3($t$-2)+3($t$-3)",
        "trans2-14": "2($t$-2)+4($t$-3)",
        "trans2-15": "1($t$-2)+5($t$-3)",
    }
    for tstr in ("m2l4", "m3l6"):
        posmap = posmaps[tstr]
        vaposmap = vaposmaps[tstr]
        mixposmap = mixposmaps[tstr]

        min_xtb_energy = min(
            i.properties["energy_xtb_dmso_au"]
            for i in atomlite.Database(
                project_dirs.directories["data_dir"] / "cage_data.db"
            ).get_entries()
            if tstr in i.key and "energy_xtb_dmso_au" in i.properties
        )

        min_xtbpbe0_energy = min(
            i.properties["energy_pbe0def2svp_dmso_au"]
            for i in atomlite.Database(
                project_dirs.directories["data_dir"] / "cage_data.db"
            ).get_entries()
            if tstr in i.key and "energy_pbe0def2svp_dmso_au" in i.properties
        )

        min_xtbwb97_energy = min(
            i.properties["energy_wb97def2tzvpp_dmso_au"]
            for i in atomlite.Database(
                project_dirs.directories["data_dir"] / "cage_data.db"
            ).get_entries()
            if tstr in i.key and "energy_wb97def2tzvpp_dmso_au" in i.properties
        )

        gxtb_entries = {
            i.key: i.properties
            for i in atomlite.Database(
                project_dirs.directories["data_dir"] / "gxtb_cage_data.db"
            ).get_entries()
        }

        min_gxtb_energy = min(
            gxtb_entries[i]["energy_gxtb_gas_au"]
            for i in gxtb_entries
            if tstr in i and "energy_gxtb_gas_au" in gxtb_entries[i]
        )

        min_gxtbpbe0_energy = min(
            gxtb_entries[i]["energy_gxtb_pbe0def2svp_gas_au"]
            for i in gxtb_entries
            if tstr in i
            and "energy_gxtb_pbe0def2svp_gas_au" in gxtb_entries[i]
        )

        min_gxtbpbe0dmso_energy = min(
            gxtb_entries[i]["energy_gxtb_pbe0def2svp_dmso_au"]
            for i in gxtb_entries
            if tstr in i
            and "energy_gxtb_pbe0def2svp_dmso_au" in gxtb_entries[i]
        )

        min_gxtbwb97_energy = min(
            gxtb_entries[i]["energy_gxtb_wb97def2tzvpp_gas_au"]
            for i in gxtb_entries
            if tstr in i
            and "energy_gxtb_wb97def2tzvpp_gas_au" in gxtb_entries[i]
        )

        # Plotting.
        fig, ax = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if tstr not in entry.key or entry.key.split("_")[-1] not in posmap:
                continue

            xtb_energy = (
                entry.properties["energy_xtb_dmso_au"] - min_xtb_energy
            ) * 2625.5
            ax.bar(
                [posmap[entry.key.split("_")[-1]] - 0.2],
                [xtb_energy],
                width=0.1,
                ec="k",
                fc="tab:blue",
                label="xTB/DMSO"
                if posmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if xtb_energy == 0:
                ax.bar(
                    [posmap[entry.key.split("_")[-1]] - 0.2],
                    [xtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = entry.properties["energy_pbe0def2svp_dmso_au"]
            except KeyError:
                continue
            pbe_energy = (
                entry.properties["energy_pbe0def2svp_dmso_au"]
                - min_xtbpbe0_energy
            ) * 2625.5
            ax.bar(
                [posmap[entry.key.split("_")[-1]]],
                [pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0/DMSO"
                if posmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if pbe_energy == 0:
                ax.bar(
                    [posmap[entry.key.split("_")[-1]]],
                    [pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = entry.properties["energy_wb97def2tzvpp_dmso_au"]
            except KeyError:
                continue
            wb97_energy = (
                entry.properties["energy_wb97def2tzvpp_dmso_au"]
                - min_xtbwb97_energy
            ) * 2625.5
            ax.bar(
                [posmap[entry.key.split("_")[-1]] + 0.2],
                [wb97_energy],
                width=0.1,
                ec="k",
                fc="tab:green",
                label="wB97x-D3/DMSO"
                if posmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if wb97_energy == 0:
                ax.bar(
                    [posmap[entry.key.split("_")[-1]] + 0.2],
                    [wb97_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax.tick_params(axis="both", which="major", labelsize=16)
        ax.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax.legend(fontsize=16)
        ax.set_xticks([posmap[i] for i in posmap])
        ax.set_xticklabels(
            [name_map[i] for i in posmap],
            fontsize=16,
            rotation=0,
        )
        ax.set_title("from GFN2-xTB", fontsize=16)
        fig.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_1_cis_v_trans.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_1_cis_v_trans.pdf",
            dpi=360,
            bbox_inches="tight",
        )

        fig2, ax1 = plt.subplots(figsize=(8, 5))

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in posmap:
                continue

            gxtb_energy = (
                gxtb_props["energy_gxtb_gas_au"] - min_gxtb_energy
            ) * 2625.5
            ax1.bar(
                [posmap[key.split("_")[-1]] - 0.3],
                [gxtb_energy],
                width=0.1,
                ec="k",
                fc="tab:cyan",
                label="g-xTB" if posmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_energy == 0:
                ax1.bar(
                    [posmap[key.split("_")[-1]] - 0.3],
                    [gxtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
                - min_gxtbpbe0_energy
            ) * 2625.5
            ax1.bar(
                [posmap[key.split("_")[-1]] - 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0" if posmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_pbe_energy == 0:
                ax1.bar(
                    [posmap[key.split("_")[-1]] - 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
                - min_gxtbpbe0dmso_energy
            ) * 2625.5
            ax1.bar(
                [posmap[key.split("_")[-1]] + 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:red",
                label="PBE0/DMSO" if posmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_pbe_energy == 0:
                ax1.bar(
                    [posmap[key.split("_")[-1]] + 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_wb97def2tzvpp_gas_au"]
            except KeyError:
                continue

            gxtb_wb97_energy = (
                gxtb_props["energy_gxtb_wb97def2tzvpp_gas_au"]
                - min_gxtbwb97_energy
            ) * 2625.5
            ax1.bar(
                [posmap[key.split("_")[-1]] + 0.3],
                [gxtb_wb97_energy],
                width=0.1,
                ec="k",
                fc="tab:green",
                label="wB97x-D3" if posmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_wb97_energy == 0:
                ax1.bar(
                    [posmap[key.split("_")[-1]] + 0.3],
                    [gxtb_wb97_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax1.tick_params(axis="both", which="major", labelsize=16)
        ax1.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax1.set_title("from g-xTB", fontsize=16)
        ax1.legend(fontsize=16)
        ax1.set_xticks([posmap[i] for i in posmap])
        ax1.set_xticklabels(
            [name_map[i] for i in posmap],
            fontsize=16,
            rotation=0,
        )
        fig2.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_2_gxtb_cis_v_trans.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig2.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_2_gxtb_cis_v_trans.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        fig3c, ax2c = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if tstr not in entry.key or entry.key.split("_")[-1] not in posmap:
                continue

            ax2c.scatter(
                [
                    posmap[entry.key.split("_")[-1]] - 0.2
                    for _ in entry.properties["geometric_analysis"][
                        "binder_binder"
                    ]
                ],
                entry.properties["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:blue",
                label="xTB" if posmap[entry.key.split("_")[-1]] == 1 else None,
            )
            ax2c.fill_between(
                [
                    posmap[entry.key.split("_")[-1]] - 0.35,
                    posmap[entry.key.split("_")[-1]] - 0.05,
                ],
                [
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                [
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                fc="tab:blue",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in posmap:
                continue

            if posmap[key.split("_")[-1]] != 0:
                ax2c.axvline(
                    x=posmap[key.split("_")[-1]] - 0.5,
                    c="gray",
                    zorder=-2,
                )

            ax2c.scatter(
                [
                    posmap[key.split("_")[-1]] + 0.2
                    for _ in gxtb_props["geometric_analysis"]["binder_binder"]
                ],
                gxtb_props["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:orange",
                label="g-xTB" if posmap[key.split("_")[-1]] == 1 else None,
            )
            ax2c.fill_between(
                [
                    posmap[key.split("_")[-1]] + 0.35,
                    posmap[key.split("_")[-1]] + 0.05,
                ],
                [
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                [
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                fc="tab:orange",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        ax2c.tick_params(axis="both", which="major", labelsize=16)
        ax2c.set_ylabel("binder-binder angle [deg]", fontsize=16)
        ax2c.legend(fontsize=16)
        ax2c.set_xticks([posmap[i] for i in posmap])
        ax2c.set_xticklabels(
            [name_map[i] for i in posmap],
            fontsize=16,
            rotation=0,
        )
        fig3c.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_3_ct_geoms.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig3c.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_3_ct_geoms.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        if tstr == "m2l4":
            continue

        fig3, ax2 = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in vaposmap
            ):
                continue

            ax2.scatter(
                [
                    vaposmap[entry.key.split("_")[-1]] - 0.2
                    for _ in entry.properties["geometric_analysis"][
                        "binder_binder"
                    ]
                ],
                entry.properties["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:blue",
                label="xTB"
                if vaposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            ax2.fill_between(
                [
                    vaposmap[entry.key.split("_")[-1]] - 0.35,
                    vaposmap[entry.key.split("_")[-1]] - 0.05,
                ],
                [
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                [
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                fc="tab:blue",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in vaposmap:
                continue

            if vaposmap[key.split("_")[-1]] != 0:
                ax2.axvline(
                    x=vaposmap[key.split("_")[-1]] - 0.5,
                    c="gray",
                    zorder=-2,
                )

            ax2.scatter(
                [
                    vaposmap[key.split("_")[-1]] + 0.2
                    for _ in gxtb_props["geometric_analysis"]["binder_binder"]
                ],
                gxtb_props["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:orange",
                label="g-xTB" if vaposmap[key.split("_")[-1]] == 1 else None,
            )
            ax2.fill_between(
                [
                    vaposmap[key.split("_")[-1]] + 0.35,
                    vaposmap[key.split("_")[-1]] + 0.05,
                ],
                [
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                [
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                fc="tab:orange",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        ax2.tick_params(axis="both", which="major", labelsize=16)
        ax2.set_ylabel("binder-binder angle [deg]", fontsize=16)
        ax2.legend(fontsize=16)
        ax2.set_xticks([vaposmap[i] for i in vaposmap])
        ax2.set_xticklabels(
            [name_map[i] for i in vaposmap],
            fontsize=16,
            rotation=90,
        )
        fig3.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_4_va_geoms.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig3.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_4_va_geoms.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        fig3b, ax2b = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in mixposmap
            ):
                continue

            ax2b.scatter(
                [
                    mixposmap[entry.key.split("_")[-1]] - 0.2
                    for _ in entry.properties["geometric_analysis"][
                        "binder_binder"
                    ]
                ],
                entry.properties["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:blue",
                label="xTB"
                if mixposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            ax2b.fill_between(
                [
                    mixposmap[entry.key.split("_")[-1]] - 0.35,
                    mixposmap[entry.key.split("_")[-1]] - 0.05,
                ],
                [
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    max(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                [
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                    min(
                        entry.properties["geometric_analysis"]["binder_binder"]
                    ),
                ],
                fc="tab:blue",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in mixposmap:
                continue

            if mixposmap[key.split("_")[-1]] != -5:  # noqa: PLR2004
                ax2b.axvline(
                    x=mixposmap[key.split("_")[-1]] - 0.5,
                    c="gray",
                    zorder=-2,
                )

            ax2b.scatter(
                [
                    mixposmap[key.split("_")[-1]] + 0.2
                    for _ in gxtb_props["geometric_analysis"]["binder_binder"]
                ],
                gxtb_props["geometric_analysis"]["binder_binder"],
                ec="k",
                fc="tab:orange",
                label="g-xTB" if mixposmap[key.split("_")[-1]] == 1 else None,
            )
            ax2b.fill_between(
                [
                    mixposmap[key.split("_")[-1]] + 0.35,
                    mixposmap[key.split("_")[-1]] + 0.05,
                ],
                [
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                    max(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                [
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                    min(gxtb_props["geometric_analysis"]["binder_binder"]),
                ],
                fc="tab:orange",
                ec="none",
                alpha=0.5,
                zorder=-2,
            )

        ax2b.tick_params(axis="both", which="major", labelsize=16)
        ax2b.set_ylabel("binder-binder angle [deg]", fontsize=16)
        ax2b.legend(fontsize=16)
        ax2b.set_xticks([mixposmap[i] for i in mixposmap])
        ax2b.set_xticklabels(
            [name_map[i] for i in mixposmap],
            fontsize=16,
            rotation=90,
        )
        fig3b.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_5_mix_geoms.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig3b.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_5_mix_geoms.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        fig4, ax4 = plt.subplots(figsize=(8, 5))

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in vaposmap:
                continue

            gxtb_energy = (
                gxtb_props["energy_gxtb_gas_au"] - min_gxtb_energy
            ) * 2625.5
            ax4.bar(
                [vaposmap[key.split("_")[-1]] - 0.3],
                [gxtb_energy],
                width=0.1,
                ec="k",
                fc="tab:cyan",
                label="g-xTB" if vaposmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_energy == 0:
                ax4.bar(
                    [vaposmap[key.split("_")[-1]] - 0.3],
                    [gxtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
                - min_gxtbpbe0_energy
            ) * 2625.5
            ax4.bar(
                [vaposmap[key.split("_")[-1]] - 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0" if vaposmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_pbe_energy == 0:
                ax4.bar(
                    [vaposmap[key.split("_")[-1]] - 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
                - min_gxtbpbe0dmso_energy
            ) * 2625.5
            ax4.bar(
                [vaposmap[key.split("_")[-1]] + 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:red",
                label="PBE0/DMSO"
                if vaposmap[key.split("_")[-1]] == 0
                else None,
            )
            if gxtb_pbe_energy == 0:
                ax4.bar(
                    [vaposmap[key.split("_")[-1]] + 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_wb97def2tzvpp_gas_au"]
            except KeyError:
                continue

            gxtb_wb97_energy = (
                gxtb_props["energy_gxtb_wb97def2tzvpp_gas_au"]
                - min_gxtbwb97_energy
            ) * 2625.5
            ax4.bar(
                [vaposmap[key.split("_")[-1]] + 0.3],
                [gxtb_wb97_energy],
                width=0.1,
                ec="k",
                fc="tab:green",
                label="wB97x-D3"
                if vaposmap[key.split("_")[-1]] == 0
                else None,
            )
            if gxtb_wb97_energy == 0:
                ax4.bar(
                    [vaposmap[key.split("_")[-1]] + 0.3],
                    [gxtb_wb97_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax4.tick_params(axis="both", which="major", labelsize=16)
        ax4.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax4.legend(fontsize=16)
        ax4.set_xticks([vaposmap[i] for i in vaposmap])
        ax4.set_xticklabels(
            [name_map[i] for i in vaposmap],
            fontsize=16,
            rotation=90,
        )
        ax4.set_title("from g-xTB", fontsize=16)
        fig4.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_6_gxtb_t2va.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig4.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_6_gxtb_t2va.pdf",
            dpi=360,
            bbox_inches="tight",
        )

        fig5, ax5 = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in vaposmap
            ):
                continue

            xtb_energy = (
                entry.properties["energy_xtb_dmso_au"] - min_xtb_energy
            ) * 2625.5
            ax5.bar(
                [vaposmap[entry.key.split("_")[-1]] - 0.1],
                [xtb_energy],
                width=0.1,
                ec="k",
                fc="tab:blue",
                label="xTB/DMSO"
                if vaposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if xtb_energy == 0:
                ax5.bar(
                    [vaposmap[entry.key.split("_")[-1]] - 0.1],
                    [xtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = entry.properties["energy_pbe0def2svp_dmso_au"]
            except KeyError:
                continue
            pbe_energy = (
                entry.properties["energy_pbe0def2svp_dmso_au"]
                - min_xtbpbe0_energy
            ) * 2625.5
            ax5.bar(
                [vaposmap[entry.key.split("_")[-1]] + 0.1],
                [pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0/DMSO"
                if vaposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if pbe_energy == 0:
                ax5.bar(
                    [vaposmap[entry.key.split("_")[-1]] + 0.1],
                    [pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax5.tick_params(axis="both", which="major", labelsize=16)
        ax5.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax5.set_title("from GFN2-xTB", fontsize=16)
        ax5.legend(fontsize=16)
        ax5.set_xticks([vaposmap[i] for i in vaposmap])
        ax5.set_xticklabels(
            [name_map[i] for i in vaposmap],
            fontsize=16,
            rotation=90,
        )
        fig5.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_7_t2va.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig5.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_7_t2va.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        fig6, ax6 = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in mixposmap
            ):
                continue

            xtb_energy = (
                entry.properties["energy_xtb_dmso_au"] - min_xtb_energy
            ) * 2625.5
            ax6.bar(
                [mixposmap[entry.key.split("_")[-1]] - 0.1],
                [xtb_energy],
                width=0.1,
                ec="k",
                fc="tab:blue",
                label="xTB/DMSO"
                if mixposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if xtb_energy == 0:
                ax6.bar(
                    [mixposmap[entry.key.split("_")[-1]] - 0.1],
                    [xtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = entry.properties["energy_pbe0def2svp_dmso_au"]
            except KeyError:
                continue
            pbe_energy = (
                entry.properties["energy_pbe0def2svp_dmso_au"]
                - min_xtbpbe0_energy
            ) * 2625.5
            ax6.bar(
                [mixposmap[entry.key.split("_")[-1]] + 0.1],
                [pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0/DMSO"
                if mixposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if pbe_energy == 0:
                ax6.bar(
                    [mixposmap[entry.key.split("_")[-1]] + 0.1],
                    [pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax6.tick_params(axis="both", which="major", labelsize=16)
        ax6.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax6.legend(fontsize=16)
        ax6.set_xticks([mixposmap[i] for i in mixposmap])
        ax6.set_xticklabels(
            [name_map[i] for i in mixposmap],
            fontsize=16,
            rotation=90,
        )
        ax6.set_title("from GFN2-xTB", fontsize=16)
        fig6.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_8_mixlig.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig6.savefig(
            project_dirs.directories["figure_dir"] / f"{tstr}_8_mixlig.pdf",
            dpi=360,
            bbox_inches="tight",
        )

        fig7, ax7 = plt.subplots(figsize=(8, 5))

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in mixposmap:
                continue

            gxtb_energy = (
                gxtb_props["energy_gxtb_gas_au"] - min_gxtb_energy
            ) * 2625.5
            ax7.bar(
                [mixposmap[key.split("_")[-1]] - 0.2],
                [gxtb_energy],
                width=0.1,
                ec="k",
                fc="tab:cyan",
                label="g-xTB" if mixposmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_energy == 0:
                ax7.bar(
                    [mixposmap[key.split("_")[-1]] - 0.2],
                    [gxtb_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_gas_au"]
                - min_gxtbpbe0_energy
            ) * 2625.5
            ax7.bar(
                [mixposmap[key.split("_")[-1]]],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="PBE0" if mixposmap[key.split("_")[-1]] == 0 else None,
            )
            if gxtb_pbe_energy == 0:
                ax7.bar(
                    [mixposmap[key.split("_")[-1]]],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
                - min_gxtbpbe0dmso_energy
            ) * 2625.5
            ax7.bar(
                [mixposmap[key.split("_")[-1]] + 0.2],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:red",
                label="PBE0/DMSO"
                if mixposmap[key.split("_")[-1]] == 0
                else None,
            )
            if gxtb_pbe_energy == 0:
                ax7.bar(
                    [mixposmap[key.split("_")[-1]] + 0.2],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax7.tick_params(axis="both", which="major", labelsize=16)
        ax7.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax7.set_title("from g-xTB", fontsize=16)
        ax7.legend(fontsize=16)
        ax7.set_xticks([mixposmap[i] for i in mixposmap])
        ax7.set_xticklabels(
            [name_map[i] for i in mixposmap],
            fontsize=16,
            rotation=90,
        )
        fig7.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_9_gxtb_mixlig.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig7.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_9_gxtb_mixlig.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        min_pbe0dmso_energy = min(
            (min_xtbpbe0_energy, min_gxtbpbe0dmso_energy)
        )

        fig8, ax8 = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in mixposmap
            ):
                continue

            try:
                _ = entry.properties["energy_pbe0def2svp_dmso_au"]
            except KeyError:
                continue
            pbe_energy = (
                entry.properties["energy_pbe0def2svp_dmso_au"]
                - min_pbe0dmso_energy
            ) * 2625.5
            ax8.bar(
                [mixposmap[entry.key.split("_")[-1]] - 0.1],
                [pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:blue",
                label="xTB|PBE0/DMSO"
                if mixposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if pbe_energy == 0:
                ax8.bar(
                    [mixposmap[entry.key.split("_")[-1]] - 0.1],
                    [pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in mixposmap:
                continue

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
                - min_pbe0dmso_energy
            ) * 2625.5
            ax8.bar(
                [mixposmap[key.split("_")[-1]] + 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="g-xTB|PBE0/DMSO"
                if mixposmap[key.split("_")[-1]] == 0
                else None,
            )
            if gxtb_pbe_energy == 0:
                ax8.bar(
                    [mixposmap[key.split("_")[-1]] + 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax8.tick_params(axis="both", which="major", labelsize=16)
        ax8.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax8.legend(fontsize=16)
        ax8.set_xticks([mixposmap[i] for i in mixposmap])
        ax8.set_xticklabels(
            [name_map[i] for i in mixposmap],
            fontsize=16,
            rotation=90,
        )
        fig8.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_10_pbe0_mixlig.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig8.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_10_pbe0_mixlig.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()

        fig9, ax9 = plt.subplots(figsize=(8, 5))

        for entry in atomlite.Database(
            project_dirs.directories["data_dir"] / "cage_data.db"
        ).get_entries():
            if (
                tstr not in entry.key
                or entry.key.split("_")[-1] not in vaposmap
            ):
                continue

            try:
                _ = entry.properties["energy_pbe0def2svp_dmso_au"]
            except KeyError:
                continue
            pbe_energy = (
                entry.properties["energy_pbe0def2svp_dmso_au"]
                - min_pbe0dmso_energy
            ) * 2625.5
            ax9.bar(
                [vaposmap[entry.key.split("_")[-1]] - 0.1],
                [pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:blue",
                label="xTB|PBE0/DMSO"
                if vaposmap[entry.key.split("_")[-1]] == 1
                else None,
            )
            if pbe_energy == 0:
                ax9.bar(
                    [vaposmap[entry.key.split("_")[-1]] - 0.1],
                    [pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        for key, gxtb_props in gxtb_entries.items():
            if tstr not in key or key.split("_")[-1] not in vaposmap:
                continue

            try:
                _ = gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
            except KeyError:
                continue

            gxtb_pbe_energy = (
                gxtb_props["energy_gxtb_pbe0def2svp_dmso_au"]
                - min_pbe0dmso_energy
            ) * 2625.5
            ax9.bar(
                [vaposmap[key.split("_")[-1]] + 0.1],
                [gxtb_pbe_energy],
                width=0.1,
                ec="k",
                fc="tab:orange",
                label="g-xTB|PBE0/DMSO"
                if vaposmap[key.split("_")[-1]] == 0
                else None,
            )
            if gxtb_pbe_energy == 0:
                ax9.bar(
                    [vaposmap[key.split("_")[-1]] + 0.1],
                    [gxtb_pbe_energy],
                    width=0.1,
                    ec="k",
                    fc="none",
                    lw=4,
                )

        ax9.tick_params(axis="both", which="major", labelsize=16)
        ax9.set_ylabel("rel. E [kJmol$^{-1}$]", fontsize=16)
        ax9.legend(fontsize=16)
        ax9.set_xticks([vaposmap[i] for i in vaposmap])
        ax9.set_xticklabels(
            [name_map[i] for i in vaposmap],
            fontsize=16,
            rotation=90,
        )
        fig9.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_11_pbe0_t2va.png",
            dpi=360,
            bbox_inches="tight",
        )
        fig9.savefig(
            project_dirs.directories["figure_dir"]
            / f"{tstr}_11_pbe0_t2va.pdf",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()


def workflow(  # noqa: C901, PLR0912, PLR0915
    cages: dict,
    project_dirs: makecages.project.Prj,
    skip_frequency: bool,  # noqa: FBT001
    solvent: str | None,
) -> None:
    """Run the cage workflow."""
    database = atomlite.Database(
        project_dirs.directories["data_dir"] / "cage_data.db"
    )
    reaction_factory = stk.DativeReactionFactory(
        stk.GenericReactionFactory(
            bond_orders={
                frozenset({stk.GenericFunctionalGroup, stk.SingleAtom}): 9,
            }
        )
    )

    with_solvent = solvent is not None

    for cage_name, cdict in cages.items():
        if database.has_entry(cage_name) and (
            f"energy_pbe0def2svp_{solvent}_au"
            in database.get_entry(cage_name).properties
        ):
            continue

        logger.info("building %s", cage_name)

        unopt_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_unopt.mol"
        )
        opt_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_optc.mol"
        )
        optxyz_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_optc.xyz"
        )
        fail_file = (
            project_dirs.directories["calculation_dir"]
            / f"{cage_name}_opt.failed"
        )

        if cdict["tstr"] in topology_codes:
            topology_code, _ = topology_codes[cdict["tstr"]]
            vertex_positions = topology_code.get_layout(
                layout_type="kamada",
                scale=10,
            )
        else:
            vertex_positions = None

        unopt_mol = stk.ConstructedMolecule(
            cdict["tfunc"](
                building_blocks=cdict["bbs"],
                # Only opt if needed.
                optimizer=stk.MCHammer(target_bond_length=2.5, num_steps=1000)
                if not opt_file.exists()
                else stk.NullOptimizer(),
                reaction_factory=reaction_factory,
                vertex_alignments=cdict["vas"],
                vertex_positions=vertex_positions,
            )
        )
        if not opt_file.exists():
            # Only write if needed.
            unopt_mol.write(unopt_file)
        unopt_mol = unopt_mol.with_structure_from_file(unopt_file)

        if not opt_file.exists() and not fail_file.exists():
            try:
                logger.info("optimising %s", cage_name)
                opt_mol = makecages.utilities.orca_cage_optimisation(
                    mol=unopt_mol,
                    name=cage_name,
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                    calc_dir=project_dirs.directories["calculation_dir"],
                    xtb_path=project_dirs.directories["xtb_path"],
                    gulp_path=project_dirs.directories["gulp_path"],
                    orca_path=project_dirs.directories["orca_path"],
                    solvent=solvent,
                    with_solvent=with_solvent,
                    metal_ff={46: "Pd4+2"},
                    skip_frequency=skip_frequency,
                )

                opt_mol.write(opt_file)
                opt_mol.write(optxyz_file)

            except (stko.ConvergenceError, stko.OptimizerError):
                logger.info("optimisation failed for %s", cage_name)
                with fail_file.open("w") as f:
                    f.write("failed")

        if fail_file.exists():
            # Add to database.
            database.update_entries(
                atomlite.Entry.from_rdkit(
                    key=cage_name,
                    molecule=unopt_mol.to_rdkit_mol(),
                    properties={"opt_failed": True, "heteroleptic": False},
                )
            )

        else:
            opt_mol = unopt_mol.with_structure_from_file(str(opt_file))

            if not skip_frequency:
                if solvent is not None:
                    freq_output_dir = (
                        project_dirs.directories["calculation_dir"]
                        / f"{cage_name}_xtbext_{solvent}"
                    )

                else:
                    freq_output_dir = (
                        project_dirs.directories["calculation_dir"]
                        / f"{cage_name}_xtbext_gas"
                    )

                has_negative_frequencies = (
                    makecages.utilities.orca_get_numerical_frequencies(
                        name=cage_name,
                        output_dir=freq_output_dir,
                    )
                )

                energy_au = makecages.utilities.read_orca_energy(
                    freq_output_dir
                )
                free_energy_au = makecages.utilities.read_orca_free_energy(
                    freq_output_dir
                )

            else:
                has_negative_frequencies = None
                if solvent is not None:
                    nofreq_output_dir = (
                        project_dirs.directories["calculation_dir"]
                        / f"{cage_name}_xtb_{solvent}"
                    )
                else:
                    nofreq_output_dir = (
                        project_dirs.directories["calculation_dir"]
                        / f"{cage_name}_xtb_gas"
                    )
                energy_au = makecages.utilities.read_xtb_energy(
                    nofreq_output_dir
                )
                free_energy_au = None

            # Get geometrical properties.
            ligands = stko.molecule_analysis.DecomposeMOC().decompose(
                molecule=opt_mol,
                metal_atom_nos=(46,),
            )
            tsa = stko.molecule_analysis.DitopicThreeSiteAnalyser()
            ligand_dict = defaultdict(list)
            for lig in ligands:
                # Define as building block with
                as_building_block = stk.BuildingBlock.init_from_molecule(
                    lig,
                    stko.functional_groups.ThreeSiteFactory(
                        smarts="[#6H1X3]~[#7X2]~[#6H1X3]"
                    ),
                )
                ligand_dict["binder_dist"].append(
                    tsa.get_binder_distance(as_building_block)
                )
                ligand_dict["binder_binder"].append(
                    tsa.get_binder_binder_angle(as_building_block)
                )
                ligand_dict["torsion"].append(
                    abs(tsa.get_binder_adjacent_torsion(as_building_block))
                )

            # Add to database.
            database.add_entries(
                atomlite.Entry.from_rdkit(
                    key=cage_name,
                    molecule=opt_mol.to_rdkit_mol(),
                    properties={
                        f"energy_xtb_{solvent}_au": energy_au,
                        f"free_energy_xtb_{solvent}_au": free_energy_au,
                        "opt_failed": False,
                        "heteroleptic": False,
                        "tstr": cdict["tstr"],
                        "charge": cdict["charge"],
                        "num_unpaired_e": cdict["num_unpaired_e"],
                        "has_negative_frequencies": has_negative_frequencies,
                        "geometric_analysis": ligand_dict,
                    },
                )
            )

            pbe0_spe_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_pbe0_{solvent}_spe"
            )
            pbe0_spe_output = pbe0_spe_output_dir / "spe.out"
            if not pbe0_spe_output.exists():
                logger.info("    ORCA PBE0/def2-SVP of %s", cage_name)

                makecages.utilities.orca_spe(
                    molecule=opt_mol,
                    output_dir=pbe0_spe_output_dir,
                    solvent=solvent,
                    method_string="PBE0 def2-SVP D3BJ",
                    pal_string="PAL8",
                    orca_path=project_dirs.directories["orca_path"],
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                )

            logger.info("loading from %s", pbe0_spe_output_dir)
            pbe_energy_au = makecages.utilities.read_orca_energy(
                pbe0_spe_output_dir,
                prefix="spe",
            )
            if pbe_energy_au is None:
                raise RuntimeError

            # Add to database.
            database.update_properties(
                atomlite.PropertyEntry(
                    key=cage_name,
                    properties={
                        f"energy_pbe0def2svp_{solvent}_au": pbe_energy_au,
                    },
                )
            )

            if "-" in cage_name and cage_name not in (
                "m3l6_trans2-1",
                "m3l6_trans2-1a",
                "m3l6_trans2-2",
            ):
                continue

            wb97_spe_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_wb97_{solvent}_spe"
            )
            logger.info("loading from %s", wb97_spe_output_dir)
            wb97_spe_output = wb97_spe_output_dir / "spe.out"
            if not wb97_spe_output.exists():
                logger.info("    ORCA wb97X-D3/def2-TZVPP of %s", cage_name)

                makecages.utilities.orca_spe(
                    molecule=opt_mol,
                    output_dir=wb97_spe_output_dir,
                    method_string="wB97X-D3 def2-TZVPP",
                    pal_string="PAL8",
                    solvent=solvent,
                    orca_path=project_dirs.directories["orca_path"],
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                )

            wb97_energy_au = makecages.utilities.read_orca_energy(
                wb97_spe_output_dir,
                prefix="spe",
            )
            if wb97_energy_au is None:
                raise RuntimeError

            # Add to database.
            database.update_properties(
                atomlite.PropertyEntry(
                    key=cage_name,
                    properties={
                        f"energy_wb97def2tzvpp_{solvent}_au": wb97_energy_au,
                    },
                )
            )


def gxtb_workflow(  # noqa: C901, PLR0912, PLR0915
    cages: dict,
    project_dirs: makecages.project.Prj,
) -> None:
    """Run the cage workflow."""
    database = atomlite.Database(
        project_dirs.directories["data_dir"] / "gxtb_cage_data.db"
    )
    reaction_factory = stk.DativeReactionFactory(
        stk.GenericReactionFactory(
            bond_orders={
                frozenset({stk.GenericFunctionalGroup, stk.SingleAtom}): 9,
            }
        )
    )

    for cage_name, cdict in cages.items():
        if database.has_entry(cage_name) and (
            "energy_gxtb_pbe0def2svp_gas_au"
            in database.get_entry(cage_name).properties
        ):
            continue

        logger.info("building %s", cage_name)

        unopt_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_unopt.mol"
        )
        opt_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_optgxtb.mol"
        )
        optxyz_file = (
            project_dirs.directories["cage_dir"] / f"{cage_name}_optgxtb.xyz"
        )
        fail_file = (
            project_dirs.directories["calculation_dir"]
            / f"{cage_name}_optgxtb.failed"
        )

        unopt_mol = stk.ConstructedMolecule(
            cdict["tfunc"](
                building_blocks=cdict["bbs"],
                # Only opt if needed.
                optimizer=stk.MCHammer(target_bond_length=2.5, num_steps=1000)
                if not opt_file.exists()
                else stk.NullOptimizer(),
                reaction_factory=reaction_factory,
                vertex_alignments=cdict["vas"],
            )
        )
        if not opt_file.exists():
            # Only write if needed.
            unopt_mol.write(unopt_file)
        unopt_mol = unopt_mol.with_structure_from_file(unopt_file)

        if not opt_file.exists() and not fail_file.exists():
            try:
                logger.info("optimising %s", cage_name)
                opt_mol = makecages.utilities.gxtb_cage_optimisation(
                    mol=unopt_mol,
                    name=cage_name,
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                    calc_dir=project_dirs.directories["calculation_dir"],
                    gxtb_path=project_dirs.directories["gxtb_path"],
                    gulp_path=project_dirs.directories["gulp_path"],
                    metal_ff={46: "Pd4+2"},
                )

                opt_mol.write(opt_file)
                opt_mol.write(optxyz_file)

            except (stko.ConvergenceError, stko.OptimizerError):
                logger.info("optimisation failed for %s", cage_name)
                with fail_file.open("w") as f:
                    f.write("failed")

        if fail_file.exists():
            # Add to database.
            database.update_entries(
                atomlite.Entry.from_rdkit(
                    key=cage_name,
                    molecule=unopt_mol.to_rdkit_mol(),
                    properties={"opt_failed": True, "heteroleptic": False},
                )
            )

        else:
            opt_mol = unopt_mol.with_structure_from_file(str(opt_file))
            nofreq_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_gxtb_gas"
            )
            energy_au = makecages.utilities.read_gxtb_energy(nofreq_output_dir)

            # Get geometrical properties.
            ligands = stko.molecule_analysis.DecomposeMOC().decompose(
                molecule=opt_mol,
                metal_atom_nos=(46,),
            )
            tsa = stko.molecule_analysis.DitopicThreeSiteAnalyser()
            ligand_dict = defaultdict(list)
            for lig in ligands:
                # Define as building block with
                as_building_block = stk.BuildingBlock.init_from_molecule(
                    lig,
                    stko.functional_groups.ThreeSiteFactory(
                        smarts="[#6H1X3]~[#7X2]~[#6H1X3]"
                    ),
                )
                ligand_dict["binder_dist"].append(
                    tsa.get_binder_distance(as_building_block)
                )
                ligand_dict["binder_binder"].append(
                    tsa.get_binder_binder_angle(as_building_block)
                )
                ligand_dict["torsion"].append(
                    abs(tsa.get_binder_adjacent_torsion(as_building_block))
                )

            # Add to database.
            database.add_entries(
                atomlite.Entry.from_rdkit(
                    key=cage_name,
                    molecule=opt_mol.to_rdkit_mol(),
                    properties={
                        "energy_gxtb_gas_au": energy_au,
                        "opt_failed": False,
                        "heteroleptic": False,
                        "tstr": cdict["tstr"],
                        "charge": cdict["charge"],
                        "num_unpaired_e": cdict["num_unpaired_e"],
                        "geometric_analysis": ligand_dict,
                    },
                )
            )

            solvent = "gas"
            pbe0_spe_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_gxtb_pbe0_{solvent}_spe"
            )
            pbe0_spe_output = pbe0_spe_output_dir / "spe.out"
            if not pbe0_spe_output.exists():
                logger.info("    ORCA PBE0/def2-SVP of %s", cage_name)

                makecages.utilities.orca_spe(
                    molecule=opt_mol,
                    output_dir=pbe0_spe_output_dir,
                    solvent=None,
                    method_string="PBE0 def2-SVP D3BJ",
                    pal_string="PAL8",
                    orca_path=project_dirs.directories["orca_path"],
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                )

            logger.info("loading from %s", pbe0_spe_output_dir)
            pbe_energy_au = makecages.utilities.read_orca_energy(
                pbe0_spe_output_dir,
                prefix="spe",
            )
            if pbe_energy_au is None:
                raise RuntimeError

            # Add to database.
            database.update_properties(
                atomlite.PropertyEntry(
                    key=cage_name,
                    properties={
                        "energy_gxtb_pbe0def2svp_gas_au": pbe_energy_au,
                    },
                )
            )

            solvent = "dmso"
            pbe0_spe_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_gxtb_pbe0_{solvent}_spe"
            )
            pbe0_spe_output = pbe0_spe_output_dir / "spe.out"
            if not pbe0_spe_output.exists():
                logger.info("    ORCA PBE0/def2-SVP of %s", cage_name)

                makecages.utilities.orca_spe(
                    molecule=opt_mol,
                    output_dir=pbe0_spe_output_dir,
                    solvent=solvent,
                    method_string="PBE0 def2-SVP D3BJ",
                    pal_string="PAL8",
                    orca_path=project_dirs.directories["orca_path"],
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                )

            logger.info("loading from %s", pbe0_spe_output_dir)
            pbe_energy_au = makecages.utilities.read_orca_energy(
                pbe0_spe_output_dir,
                prefix="spe",
            )
            if pbe_energy_au is None:
                raise RuntimeError

            # Add to database.
            database.update_properties(
                atomlite.PropertyEntry(
                    key=cage_name,
                    properties={
                        f"energy_gxtb_pbe0def2svp_{solvent}_au": pbe_energy_au,
                    },
                )
            )

            if "-" in cage_name and cage_name not in (
                "m3l6_trans2-1",
                "m3l6_trans2-1a",
                "m3l6_trans2-2",
            ):
                continue

            solvent = "gas"
            wb97_spe_output_dir = (
                project_dirs.directories["calculation_dir"]
                / f"{cage_name}_gxtb_wb97_{solvent}_spe"
            )
            logger.info("loading from %s", wb97_spe_output_dir)
            wb97_spe_output = wb97_spe_output_dir / "spe.out"
            if not wb97_spe_output.exists():
                logger.info("    ORCA wb97X-D3/def2-TZVPP of %s", cage_name)

                makecages.utilities.orca_spe(
                    molecule=opt_mol,
                    output_dir=wb97_spe_output_dir,
                    solvent=None,
                    method_string="wB97X-D3 def2-TZVPP",
                    pal_string="PAL8",
                    orca_path=project_dirs.directories["orca_path"],
                    charge=cdict["charge"],
                    num_unpaired_e=cdict["num_unpaired_e"],
                )

            wb97_energy_au = makecages.utilities.read_orca_energy(
                wb97_spe_output_dir,
                prefix="spe",
            )
            if wb97_energy_au is None:
                raise RuntimeError

            # Add to database.
            database.update_properties(
                atomlite.PropertyEntry(
                    key=cage_name,
                    properties={
                        "energy_gxtb_wb97def2tzvpp_gas_au": wb97_energy_au,
                    },
                )
            )


def main() -> None:  # noqa: D103
    project_dirs = makecages.project.Prj.init_from_json(
        pathlib.Path(__file__).resolve().parent / "project.json"
    )

    conformer_dir = project_dirs.directories["calculation_dir"] / "l1_confs"

    trans_building_block = stk.BuildingBlock.init_from_file(
        conformer_dir / "c185.mol",
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6H1X3]~[#7X2]~[#6H1X3]",
            bonders=(1,),
            deleters=(),
        ),
    )

    trans_building_block_2 = stk.BuildingBlock.init_from_file(
        conformer_dir / "c0_custom.mol",
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6H1X3]~[#7X2]~[#6H1X3]",
            bonders=(1,),
            deleters=(),
        ),
    )

    trans_building_block_3 = stk.BuildingBlock.init_from_file(
        conformer_dir / "c0_custom_2.mol",
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6H1X3]~[#7X2]~[#6H1X3]",
            bonders=(1,),
            deleters=(),
        ),
    )

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

    cages = {
        "m2l4_cis": {
            "bbs": {
                _pd_bb: (0, 1),
                cis_building_block: (2, 3, 4, 5),
            },
            "charge": 4,
            "num_unpaired_e": 0,
            "tstr": "m2l4",
            "tfunc": stk.cage.M2L4Lantern,
            "vas": None,
        },
        "m2l4_trans": {
            "bbs": {
                _pd_bb: (0, 1),
                trans_building_block: (2, 3, 4, 5),
            },
            "charge": 4,
            "num_unpaired_e": 0,
            "tstr": "m2l4",
            "tfunc": stk.cage.M2L4Lantern,
            "vas": None,
        },
        "m3l6_cis": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                cis_building_block: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans3": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_3: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-1": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 1, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0},
        },
        "m3l6_trans2-1a": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6a",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 1, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0},
        },
        "m3l6_trans2-2": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 0, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0},
        },
        "m3l6_trans2-3": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 0, 4: 0, 5: 0, 6: 0, 7: 1, 8: 0},
        },
        "m3l6_trans2-4": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (4, 5, 6, 7, 8),
                trans_building_block: (3,),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-5": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (4, 5, 6, 7, 8),
                trans_building_block_3: (3,),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-6": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (5, 6, 7, 8),
                trans_building_block: (3, 4),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-7": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 1, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0},
        },
        "m3l6_trans2-8": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (3, 4, 5, 6, 7, 8),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": {3: 1, 4: 1, 5: 1, 6: 0, 7: 0, 8: 0},
        },
        "m3l6_trans2-9": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (6, 7, 8),
                trans_building_block: (3, 4, 5),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-10": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (7, 8),
                trans_building_block: (3, 4, 5, 6),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-11": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (8,),
                trans_building_block: (3, 4, 5, 6, 7),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-12": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (5, 6, 7, 8),
                trans_building_block_3: (3, 4),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-13": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (6, 7, 8),
                trans_building_block_3: (3, 4, 5),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-14": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (7, 8),
                trans_building_block_3: (3, 4, 5, 6),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
        "m3l6_trans2-15": {
            "bbs": {
                _pd_bb: (0, 1, 2),
                trans_building_block_2: (8,),
                trans_building_block_3: (3, 4, 5, 6, 7),
            },
            "charge": 6,
            "num_unpaired_e": 0,
            "tstr": "m3l6",
            "tfunc": stk.cage.M3L6,
            "vas": None,
        },
    }
    workflow(
        cages=cages,
        project_dirs=project_dirs,
        solvent="dmso",
        skip_frequency=False,
    )
    gxtb_workflow(cages=cages, project_dirs=project_dirs)
    plotting(project_dirs)


if __name__ == "__main__":
    main()
