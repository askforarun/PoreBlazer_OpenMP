#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


TWO_LETTER = {
    "BR", "CL", "NA", "MG", "AL", "SI", "CA", "SC", "TI", "CR", "MN",
    "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "SR", "ZR",
    "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB",
    "TE", "CS", "BA", "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "RE", "OS", "IR",
    "PT", "AU", "HG", "TL", "PB", "BI",
}


def infer_element(atom_name: str) -> str:
    letters = "".join(ch for ch in atom_name if ch.isalpha()).upper()
    if not letters:
        raise ValueError(f"Could not infer element from atom name {atom_name!r}")
    if len(letters) >= 2 and letters[:2] in TWO_LETTER:
        return letters[:2].capitalize()
    return letters[0]


def convert_gro_to_xyz(gro_path: Path, xyz_path: Path) -> None:
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{gro_path} is too short to be a valid .gro file")

    natoms = int(lines[1].strip())
    atom_lines = lines[2:2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError(f"{gro_path} declares {natoms} atoms but has {len(atom_lines)} atom records")

    out_lines = [str(natoms), "Converted from GROMACS .gro for PoreBlazer"]
    for line in atom_lines:
        atom_name = line[10:15].strip()
        x = float(line[20:28]) * 10.0
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        element = infer_element(atom_name)
        out_lines.append(f"{element:<2s} {x:16.9f} {y:16.9f} {z:16.9f}")

    xyz_path.write_text("\n".join(out_lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a GROMACS .gro structure to a PoreBlazer-style .xyz file.")
    parser.add_argument("gro", type=Path, help="Input .gro file")
    parser.add_argument("xyz", type=Path, help="Output .xyz file")
    args = parser.parse_args()
    convert_gro_to_xyz(args.gro, args.xyz)


if __name__ == "__main__":
    main()
