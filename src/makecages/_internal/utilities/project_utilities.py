"""Utilities for project definition."""

import json
import pathlib
from dataclasses import dataclass


@dataclass(frozen=True, kw_only=True)
class Prj:
    """Define the project files and directories."""

    wd: pathlib.Path
    directories: dict[str, pathlib.Path]

    @classmethod
    def init_from_json(cls, json_path: pathlib.Path) -> "Prj":
        """Load Prj from a json."""
        with json_path.open("r") as f:
            data = json.load(f)

        other_directories = {
            i: j for i, j in data.items() if i not in ("wd", "under_wd")
        }
        directories = {}
        for name, path_str in other_directories.items():
            if not isinstance(path_str, str):
                raise TypeError

            directories[name] = pathlib.Path(path_str)

        for name, path_list in data["under_wd"].items():
            if not isinstance(path_list, list):
                raise TypeError

            path = pathlib.Path(data["wd"]) / "/".join(path_list)
            directories[name] = path
            path.mkdir(exist_ok=True, parents=True)

        return Prj(wd=pathlib.Path(data["wd"]), directories=directories)
