"""Utilities for Orca calculations."""

import logging
import os
import pathlib
import shutil
import subprocess as sp

import stk

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def run_job(
    molecule: stk.Molecule,
    output_path: pathlib.Path,
    orca_path: pathlib.Path,
    input_file_string: str,
    prefix: str,
) -> None:
    """Run a calculationin Orca."""
    if output_path.exists():
        shutil.rmtree(output_path)

    output_path.mkdir(parents=True)
    init_dir = pathlib.Path.cwd()
    os.chdir(output_path)

    try:
        # Bring the file in.
        molecule.write("initial.xyz")

        # Write the .inp.
        inp_file = pathlib.Path(f"{prefix}.inp")
        with inp_file.open("w") as f:
            f.write(input_file_string)
        # Run.
        output_file = pathlib.Path(f"{prefix}.out")
        err_file = pathlib.Path(f"{prefix}.err")
        cmd = [orca_path, inp_file]
        with output_file.open("w") as f_out, err_file.open("w") as f_err:
            _ = sp.run(  # noqa: S603
                cmd,
                check=False,
                input=None,
                stdout=f_out,
                stderr=f_err,
                text=True,
            )

        if prefix == "numfreq":
            hess_file = pathlib.Path(f"{prefix}.hess")
            if not hess_file.exists():
                msg = "Hess file not found."
                raise RuntimeError(msg)

        with output_file.open("r") as f:
            lastline = list(f)[-1]
        time_taken = lastline.strip()

        logger.info("job took %s", time_taken)

    finally:
        os.chdir(init_dir)


def orca_get_numerical_frequencies(
    name: str,
    output_dir: pathlib.Path,
) -> bool:
    """Run XTB2 numerical frequency calculation."""
    logger.info("load from\n %s", output_dir)

    output_file = output_dir / "numfreq.out"
    hess_file = output_dir / "numfreq.hess"

    if not output_file.exists():
        msg = "ORCA numerical frequencies not run."
        raise FileNotFoundError(msg)

    # Extract frequencies from .hess.
    frequencies = []
    with hess_file.open("r") as f:
        collect = False
        for line in f.readlines():
            if collect:
                frequencies.append(line.strip())
                # Check for final one.
                if (
                    len(frequencies) > 1
                    and int(line.strip().split()[0]) == int(frequencies[0]) - 1
                ):
                    break
            if "$vibrational_frequencies" in line.strip():
                collect = True

    cleaned_frequencies = [float(i.split()[1]) for i in frequencies[1:]]
    negative_frequencies = [i for i in cleaned_frequencies if i < 0]
    if len(negative_frequencies) > 0:
        logger.info("negative frequencies found for %s", name)
        neg_freq_file = output_dir.parent / f"{name}_orca_numfreq.negative"
        with neg_freq_file.open("w") as f:
            f.write(
                f"{len(negative_frequencies)} of {len(cleaned_frequencies)}"
            )
        has_negative_frequencies = True

    else:
        has_negative_frequencies = False

    return has_negative_frequencies


def orca_numerical_frequencies(  # noqa: PLR0913
    molecule: stk.Molecule,
    opt: str | None,
    output_dir: pathlib.Path,
    orca_path: pathlib.Path,
    xtb_path: pathlib.Path,
    charge: int,
    num_unpaired_e: int,
    solvent: str | None,
    temperature: float,
) -> stk.Molecule | None:
    """Run XTB2 numerical frequency calculation."""
    if not orca_path.exists():
        msg = "Orca path is not right!"
        raise FileNotFoundError(msg)

    if not xtb_path.exists():
        msg = "xtb path is not right!"
        raise FileNotFoundError(msg)

    os.environ["XTBEXE"] = str(xtb_path)

    logger.info("saving all calculations into\n %s", output_dir)

    # Settings.
    # ORCA options.
    solvent_string = "" if solvent is None else f"ALPB({solvent})"

    spin = num_unpaired_e / 2
    multiplicity = int(2 * spin + 1)

    run_type = f"{opt} NumFreq" if opt is not None else "NumFreq"

    input_file_string = (
        f"! XTB2 {run_type} TightSCF PAL4 {solvent_string}\n\n"
        f"%FREQ TEMP {temperature} END\n\n"
        f"* xyzfile {charge} {multiplicity} initial.xyz\n"
    )

    if run_type not in ("VeryTightOpt NumFreq", "NumFreq"):
        raise NotImplementedError

    if ("NumFreq" in run_type or "Freq" in run_type) and temperature is None:
        msg = "Set 'temp' to compute thermochemistry."
        raise ValueError(msg)

    output_file = output_dir / "numfreq.out"
    hess_file = output_dir / "numfreq.hess"
    xyz_file = output_dir / "numfreq.xyz"

    if not output_file.exists() or not hess_file.exists():
        logger.info(
            "running ORCA with input file string:\n %s",
            input_file_string,
        )
        run_job(
            molecule=molecule,
            output_path=output_dir,
            orca_path=orca_path,
            input_file_string=input_file_string,
            prefix="numfreq",
        )
    if opt is not None:
        return molecule.with_structure_from_file(xyz_file)
    return None


def orca_spe(  # noqa: PLR0913
    molecule: stk.Molecule,
    output_dir: pathlib.Path,
    orca_path: pathlib.Path,
    pal_string: str,
    charge: int,
    num_unpaired_e: int,
    solvent: str | None,
    method_string: str = "PBE0 def2-SVP D3BJ",
) -> None:
    """Run PBE0 SPE calculation."""
    if not orca_path.exists():
        msg = "Orca path is not right!"
        raise FileNotFoundError(msg)

    logger.info("saving all calculations into\n %s", output_dir)

    # Settings.
    # ORCA options.
    solvent_string = "" if solvent is None else f"CPCM({solvent})"

    spin = num_unpaired_e / 2
    multiplicity = int(2 * spin + 1)

    input_file_string = (
        f"! {method_string} SP TightSCF {pal_string} {solvent_string}\n\n"
        f"* xyzfile {charge} {multiplicity} initial.xyz\n"
    )

    output_file = output_dir / "spe.out"

    if not output_file.exists():
        logger.info(
            "running ORCA with input file string:\n %s",
            input_file_string,
        )
        run_job(
            molecule=molecule,
            output_path=output_dir,
            orca_path=orca_path,
            input_file_string=input_file_string,
            prefix="spe",
        )
