"""Utilities for cage processing."""

import logging
import pathlib

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)


def read_xtb_energy(output_dir: pathlib.Path) -> float | None:
    """Read xTB energy from optimisation output."""
    with (output_dir / "optimization_1.output").open("r") as f:
        for line in f.readlines():
            if "| TOTAL ENERGY" in line:
                return float(line.strip().split()[3])
    return None


def read_gxtb_energy(output_dir: pathlib.Path) -> float | None:
    """Read xTB energy from optimisation output."""
    with (output_dir / "optimization_1.output").open("r") as f:
        for line in f.readlines():
            if "| TOTAL ENERGY" in line:
                return float(line.strip().split()[3])

    return None


def read_orca_energy(
    output_dir: pathlib.Path,
    prefix: str = "numfreq",
) -> float | None:
    """Read final total energy from ORCA output."""
    with (output_dir / f"{prefix}.out").open("r") as f:
        energies = [
            float(line.strip().split()[-1])
            for line in f.readlines()
            if "FINAL SINGLE POINT ENERGY" in line
        ]

    if len(energies) > 0:
        return energies[-1]

    return None


def read_orca_free_energy(output_dir: pathlib.Path) -> float | None:
    """Read final total energy from ORCA output."""
    with (output_dir / "numfreq.out").open("r") as f:
        for line in f.readlines():
            if "Final Gibbs free energy" in line:
                return float(line.strip().split()[-2])

    return None
