"""Utilities for cage processing."""

import logging
import os
import pathlib
import subprocess as sp

import stk
import stko

from .orca_utilities import orca_numerical_frequencies

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


class gXTBEnergy(stko.XTBEnergy):
    """Uses g-xtb to calculate energy and other properties."""

    def __init__(  # noqa: PLR0913
        self,
        xtb_path: pathlib.Path | str,
        output_dir: pathlib.Path | str | None = None,
        num_cores: int = 1,
        calculate_free_energy: bool = False,
        calculate_ip_and_ea: bool = False,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        solvent_grid: str = "normal",
        charge: int = 0,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,
        write_sasa_info: bool = False,
    ) -> None:
        if solvent is not None:
            msg = "solvent not implemented for g-xtb"
            raise NotImplementedError(msg)

        xtb_path = pathlib.Path(xtb_path)
        self._check_path(xtb_path)
        self._xtb_path = xtb_path
        self._output_dir = (
            None if output_dir is None else pathlib.Path(output_dir)
        )
        self._num_cores = str(num_cores)
        self._calculate_free_energy = calculate_free_energy
        self._calculate_ip_and_ea = calculate_ip_and_ea
        self._electronic_temperature = str(electronic_temperature)
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self._solvent = solvent
        self._solvent_model = solvent_model
        self._solvent_grid = solvent_grid
        self._write_sasa_info = write_sasa_info

    def _run_xtb(
        self,
        xyz: pathlib.Path,
        out_file: pathlib.Path,
        init_dir: pathlib.Path,
        output_dir: pathlib.Path,
    ) -> None:
        """Runs GFN-xTB.

        Parameters:
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

            init_dir:
                The name of the current working directory.

            output_dir:
                The name of the directory into which files generated during
                the calculation are written.

        """
        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        if self._solvent is not None:
            raise NotImplementedError
            solvent = f"--{self._solvent_model} {self._solvent} "
        solvent = ""

        hess = "--hess" if self._calculate_free_energy else ""

        vipea = "--vipea" if self._calculate_ip_and_ea else ""

        cmd = (
            f"{memory} {self._xtb_path} "
            f"{xyz} "
            f"--gxtb "
            f"{hess} {vipea} --parallel {self._num_cores} "
            f"--etemp {self._electronic_temperature} "
            f"{solvent} --chrg {self._charge} "
            f"--uhf {self._num_unpaired_electrons} -I det_control.in"
        )

        try:
            os.chdir(output_dir)
            self._write_detailed_control()
            with out_file.open("w") as f:
                # Note that sp.call will hold the program until
                # completion of the calculation.
                sp.call(  # noqa: S602
                    cmd,
                    stdin=sp.PIPE,
                    stdout=f,
                    stderr=sp.PIPE,
                    # Shell is required to run complex arguments.
                    shell=True,
                )
        finally:
            os.chdir(init_dir)


class gXTB(stko.XTB):
    """Uses g-xTB to optimize molecules."""

    def __init__(  # noqa: PLR0913
        self,
        xtb_path: pathlib.Path | str,
        output_dir: pathlib.Path | str | None = None,
        opt_level: str = "normal",
        max_runs: int = 2,
        calculate_hessian: bool = True,
        num_cores: int = 1,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        solvent_grid: str = "normal",
        charge: int = 0,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,
        write_sasa_info: bool = False,
    ) -> None:
        if solvent is not None:
            msg = "solvent not implemented for g-xtb"
            raise NotImplementedError(msg)

        if not calculate_hessian and max_runs != 1:
            max_runs = 1
            logger.warning(
                "Requested that hessian calculation was skipped "
                "but the number of optimizations requested was "
                "greater than 1. The number of optimizations has been "
                "set to 1."
            )

        self._check_path(xtb_path)
        self._xtb_path = pathlib.Path(xtb_path)
        self._output_dir = (
            None if output_dir is None else pathlib.Path(output_dir)
        )
        self._opt_level = opt_level
        self._max_runs = max_runs
        self._calculate_hessian = calculate_hessian
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_model = solvent_model
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self._write_sasa_info = write_sasa_info
        self.incomplete = set()

    def _run_xtb(self, xyz: str, out_file: pathlib.Path | str) -> None:
        """Run GFN-xTB.

        Parameters
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

        """
        out_file = pathlib.Path(out_file)

        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        # Set optimization level and type.
        if self._calculate_hessian:
            # Do optimization and check hessian.
            optimization = f"--ohess {self._opt_level}"
        else:
            # Do optimization.
            optimization = f"--opt {self._opt_level}"

        if self._solvent is not None:
            raise NotImplementedError
            solvent = f"--{self._solvent_model} {self._solvent} "
        solvent = ""

        cmd = (
            f"{memory} {self._xtb_path} {xyz} "
            f"--gxtb "
            f"{optimization} --parallel {self._num_cores} "
            f"--etemp {self._electronic_temperature} "
            f"{solvent} --chrg {self._charge} "
            f"--uhf {self._num_unpaired_electrons} -I det_control.in"
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )


def gxtb_cage_optimisation(  # noqa: PLR0913
    mol: stk.Molecule,
    name: str,
    charge: int,
    num_unpaired_e: int,
    calc_dir: pathlib.Path,
    gxtb_path: pathlib.Path,
    gulp_path: pathlib.Path,
    metal_ff: dict[int, str],
) -> stk.Molecule:
    """Run optimisation workflow."""
    gulp1_output = calc_dir / f"{name}_gulp1.mol"
    gulp2_output = calc_dir / f"{name}_gulp2.mol"
    gxtb1_output = calc_dir / f"{name}_gxtbnomd.mol"

    if not gxtb_path.exists():
        msg = f"xtb is not installed here: {gxtb_path}"
        raise ValueError(msg)
    if not gulp_path.exists():
        msg = f"gulp is not installed here: {gulp_path}"
        raise ValueError(msg)

    if not gulp1_output.exists():
        output_dir = calc_dir / f"{name}_gulp1"

        logger.info("    UFF4MOF optimisation 1 of %s CG: True", name)
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path,
            maxcyc=1000,
            metal_FF=metal_ff,
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=True,
        )
        gulp_opt.assign_FF(mol)
        gulp1_mol = gulp_opt.optimize(mol=mol)
        gulp1_mol.write(gulp1_output)
    else:
        logger.info("    loading %s", gulp1_output)
        gulp1_mol = mol.with_structure_from_file(str(gulp1_output))

    if not gulp2_output.exists():
        output_dir = calc_dir / f"{name}_gulp2"

        logger.info("    UFF4MOF optimisation 2 of %s CG: False", name)
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path,
            maxcyc=1000,
            metal_FF=metal_ff,
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=False,
        )
        gulp_opt.assign_FF(gulp1_mol)
        gulp2_mol = gulp_opt.optimize(mol=gulp1_mol)
        gulp2_mol.write(gulp2_output)
    else:
        logger.info("    loading %s", gulp2_output)
        gulp2_mol = mol.with_structure_from_file(str(gulp2_output))

    if not gxtb1_output.exists():
        output_dir = calc_dir / f"{name}_gxtb_gas"

        logger.info("    g-xtb opt of %s", name)
        xtb_opt = gXTB(
            xtb_path=gxtb_path,
            output_dir=output_dir,
            num_cores=6,
            charge=charge,
            opt_level="normal",
            num_unpaired_electrons=num_unpaired_e,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True,
            solvent=None,
        )
        gxtb1_mol = xtb_opt.optimize(mol=gulp2_mol)
        gxtb1_mol.write(gxtb1_output)

    else:
        logger.info("    loading %s", gxtb1_output)
        gxtb1_mol = mol.with_structure_from_file(str(gxtb1_output))

    return mol.with_structure_from_file(str(gxtb1_output))


def orca_cage_optimisation(  # noqa: C901, PLR0912, PLR0913, PLR0915
    mol: stk.Molecule,
    name: str,
    charge: int,
    num_unpaired_e: int,
    calc_dir: pathlib.Path,
    xtb_path: pathlib.Path,
    gulp_path: pathlib.Path,
    orca_path: pathlib.Path,
    metal_ff: dict[int, str],
    solvent: str | None,
    skip_frequency: bool,  # noqa: FBT001
    with_solvent: bool = True,  # noqa: FBT001, FBT002
) -> stk.Molecule:
    """Run optimisation workflow."""
    gulp1_output = calc_dir / f"{name}_gulp1.mol"
    gulp2_output = calc_dir / f"{name}_gulp2.mol"
    if solvent is not None:
        xtbopt2_output = calc_dir / f"{name}_xtbnomd_{solvent}.mol"
        xtbopt3_output = calc_dir / f"{name}_xtbextnomd_{solvent}.mol"

    else:
        xtbopt2_output = calc_dir / f"{name}_xtb2.mol"
        xtbopt3_output = calc_dir / f"{name}_xtb3.mol"

    if not xtb_path.exists():
        msg = f"xtb is not installed here: {xtb_path}"
        raise ValueError(msg)
    if not gulp_path.exists():
        msg = f"gulp is not installed here: {gulp_path}"
        raise ValueError(msg)

    if solvent is None and with_solvent:
        msg = "Solvent must be defined."
        raise ValueError(msg)

    if solvent is not None and not with_solvent:
        msg = "Solvent defined, but not requested."
        raise ValueError(msg)

    if not gulp1_output.exists():
        output_dir = calc_dir / f"{name}_gulp1"

        logger.info("    UFF4MOF optimisation 1 of %s CG: True", name)
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path,
            maxcyc=1000,
            metal_FF=metal_ff,
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=True,
        )
        gulp_opt.assign_FF(mol)
        gulp1_mol = gulp_opt.optimize(mol=mol)
        gulp1_mol.write(gulp1_output)
    else:
        logger.info("    loading %s", gulp1_output)
        gulp1_mol = mol.with_structure_from_file(str(gulp1_output))

    if not gulp2_output.exists():
        output_dir = calc_dir / f"{name}_gulp2"

        logger.info("    UFF4MOF optimisation 2 of %s CG: False", name)
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path,
            maxcyc=1000,
            metal_FF=metal_ff,
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=False,
        )
        gulp_opt.assign_FF(gulp1_mol)
        gulp2_mol = gulp_opt.optimize(mol=gulp1_mol)
        gulp2_mol.write(gulp2_output)
    else:
        logger.info("    loading %s", gulp2_output)
        gulp2_mol = mol.with_structure_from_file(str(gulp2_output))

    if skip_frequency:
        raise NotImplementedError

    if not xtbopt3_output.exists():
        if with_solvent:
            output_dir = calc_dir / f"{name}_xtbext_{solvent}"
        else:
            output_dir = calc_dir / f"{name}_xtbext_gas"

        logger.info("    ORCA VeryTightOpt|NumFreq of %s", name)

        xtbextsolvopt_mol = orca_numerical_frequencies(
            molecule=gulp2_mol,
            opt="VeryTightOpt",
            output_dir=output_dir,
            solvent=solvent,
            orca_path=orca_path,
            xtb_path=xtb_path,
            charge=charge,
            temperature=293.15,  # K
            num_unpaired_e=num_unpaired_e,
        )
        xtbextsolvopt_mol.write(xtbopt3_output)

    else:
        logger.info("    loading %s", xtbopt3_output)
        xtbextsolvopt_mol = mol.with_structure_from_file(str(xtbopt3_output))

    return mol.with_structure_from_file(str(xtbopt3_output))
