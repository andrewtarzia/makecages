"""Build ligands."""

import json
import logging
import pathlib

import bbprep
import stk
import stko

import makecages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def xtb_geometry_energy(  # noqa: PLR0913
    building_block: stk.BuildingBlock,
    xtb_path: pathlib.Path,
    calculation_dir: pathlib.Path,
    molecule_name: str,
    charge: int,
    num_unpaired_e: int,
    save_dir: pathlib.Path,
) -> stk.BuildingBlock:
    """Geometry optimise and get energy with GFN2-xTB."""
    building_block = stko.XTB(
        xtb_path=xtb_path,
        output_dir=calculation_dir / f"{molecule_name}",
        max_runs=1,
        calculate_hessian=False,
        charge=charge,
        num_unpaired_electrons=num_unpaired_e,
        num_cores=4,
        solvent_model="alpb",
        solvent="dmso",
        solvent_grid="verytight",
        unlimited_memory=True,
    ).optimize(building_block)

    energy = stko.XTBEnergy(
        xtb_path=xtb_path,
        output_dir=calculation_dir / f"{molecule_name}_ey",
        num_cores=4,
        charge=charge,
        unlimited_memory=True,
        num_unpaired_electrons=num_unpaired_e,
        solvent_model="alpb",
        solvent="dmso",
        solvent_grid="verytight",
    ).get_energy(building_block)

    with (save_dir / f"{molecule_name}_optl.ey").open("w") as f:
        f.write(str(energy))

    return building_block


def load_ligands(json_file: pathlib.Path) -> dict[str, dict]:
    """Load ligands from json."""
    with json_file.open("r") as f:
        return json.load(f)


def main() -> None:  # noqa: D103, PLR0915
    project_dirs = makecages.project.Prj.init_from_json(
        pathlib.Path(__file__).resolve().parent / "project.json"
    )

    molecules = load_ligands(
        project_dirs.directories["ligand_dir"] / "molecules.json"
    )

    for molecule_name, mdict in molecules.items():
        if mdict["type"] in ("ion", "solvent"):
            save_dir = project_dirs.directories["ion_dir"]
            do_conformer_search = False
        else:
            save_dir = project_dirs.directories["ligand_dir"]
            do_conformer_search = True

        opt_file = save_dir / f"{molecule_name}_optl.mol"

        if not opt_file.exists():
            logger.info("optimising %s", molecule_name)
            if not do_conformer_search:
                building_block = stk.BuildingBlock(mdict["smiles"])
                building_block = xtb_geometry_energy(
                    building_block=building_block,
                    xtb_path=project_dirs.directories["xtb_path"],
                    calculation_dir=project_dirs.directories[
                        "calculation_dir"
                    ],
                    molecule_name=molecule_name,
                    num_unpaired_e=0,
                    charge=mdict["charge"],
                    save_dir=save_dir,
                )
                building_block.write(save_dir / f"{molecule_name}_optl.mol")

                energy = stko.XTBEnergy(
                    xtb_path=project_dirs.directories["xtb_path"],
                    output_dir=project_dirs.directories["calculation_dir"]
                    / f"{molecule_name}_ey",
                    num_cores=4,
                    charge=mdict["charge"],
                    unlimited_memory=True,
                    solvent_model="alpb",
                    solvent="dmso",
                    solvent_grid="verytight",
                ).get_energy(building_block)

                with (save_dir / f"{molecule_name}_optl.ey").open("w") as f:
                    f.write(str(energy))

            else:
                building_block = stk.BuildingBlock(
                    smiles=mdict["smiles"],
                    functional_groups=(
                        stk.SmartsFunctionalGroupFactory(
                            smarts="[#6]~[#7X2]~[#6]",
                            bonders=(1,),
                            deleters=(),
                        ),
                    ),
                )

                generator = bbprep.generators.ETKDG(num_confs=200)
                ensemble = generator.generate_conformers(building_block)
                ensemble = ensemble.optimise_conformers(
                    optimiser=bbprep.Optimiser(
                        name="MMFF",
                        function=stko.MMFF().optimize,
                    ),
                )
                ensemble = ensemble.optimise_conformers(
                    optimiser=bbprep.Optimiser(
                        name="xtb",
                        function=stko.XTB(
                            xtb_path=project_dirs.directories["xtb_path"],
                            output_dir=project_dirs.directories[
                                "calculation_dir"
                            ]
                            / f"{molecule_name}",
                            max_runs=1,
                            calculate_hessian=False,
                            num_cores=4,
                            unlimited_memory=True,
                        ).optimize,
                    ),
                )

                (
                    project_dirs.directories["calculation_dir"]
                    / f"{molecule_name}_confs"
                ).mkdir(exist_ok=True)
                for conformer in ensemble.yield_conformers():
                    conformer.molecule.write(
                        project_dirs.directories["calculation_dir"]
                        / f"{molecule_name}_confs"
                        / f"c{conformer.conformer_id}.mol"
                    )

                minimum_conformer = ensemble.get_lowest_energy_conformer(
                    calculator=bbprep.EnergyCalculator(
                        name="xtbenergy",
                        function=stko.XTBEnergy(
                            xtb_path=project_dirs.directories["xtb_path"],
                            output_dir=project_dirs.directories[
                                "calculation_dir"
                            ]
                            / f"{molecule_name}_ey",
                            num_cores=4,
                            unlimited_memory=True,
                        ).get_energy,
                    ),
                )
                minimum_conformer.molecule.write(
                    save_dir / f"{molecule_name}_optl.mol"
                )

        chosen = {
            "c185": "t1",
            "c0_custom": "t2",
            "c0_custom_2": "t3",
            "c63": "cis",
        }

        energies = {}
        gxtb_energies = {}
        for molfile in (
            project_dirs.directories["calculation_dir"]
            / f"{molecule_name}_confs"
        ).glob("*.mol"):
            confname = molfile.stem
            confmol = stk.BuildingBlock.init_from_file(molfile)

            eyfile = (
                project_dirs.directories["calculation_dir"]
                / f"{confname}.xtbey"
            )
            if not eyfile.exists():
                ey = stko.XTBEnergy(
                    xtb_path=project_dirs.directories["xtb_path"],
                    output_dir=project_dirs.directories["calculation_dir"]
                    / f"{confname}_ey",
                    num_cores=4,
                    unlimited_memory=True,
                ).get_energy(confmol)
                with eyfile.open("w") as f:
                    f.write(str(ey))

            with eyfile.open("r") as f:
                ey = float(f.readlines()[0].strip())
            energies[confname] = ey

            eyfile = (
                project_dirs.directories["calculation_dir"]
                / f"{confname}.gxtbey"
            )
            if not eyfile.exists():
                ey = makecages.utilities.gXTBEnergy(
                    xtb_path=project_dirs.directories["gxtb_path"],
                    output_dir=project_dirs.directories["calculation_dir"]
                    / f"{confname}_gxtb_ey",
                    num_cores=4,
                    unlimited_memory=True,
                ).get_energy(confmol)
                with eyfile.open("w") as f:
                    f.write(str(ey))

            with eyfile.open("r") as f:
                ey = float(f.readlines()[0].strip())
            gxtb_energies[confname] = ey

        energies = {
            i: (j - min(energies.values())) * 2625.5
            for i, j in energies.items()
        }
        gxtb_energies = {
            i: (j - min(gxtb_energies.values())) * 2625.5
            for i, j in gxtb_energies.items()
        }

        logger.info(
            "minimums: %s and %s (g-xtb)",
            next(i for i, j in energies.items() if j == 0),
            next(i for i, j in gxtb_energies.items() if j == 0),
        )
        for i, purpose in chosen.items():
            logger.info(
                "chosen %s for %s: %s and %s (g-xtb)",
                i,
                purpose,
                round(energies[i], 2),
                round(gxtb_energies[i], 2),
            )


if __name__ == "__main__":
    main()
