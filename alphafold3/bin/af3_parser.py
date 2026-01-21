# BioFeatureFactory
# Copyright (C) 2023–2026  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
AlphaFold3 output parser.

Parses mmCIF structure files and JSON confidence files to extract
binding metrics for RNA-protein complexes.
"""

import json
import sys
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import math


@dataclass
class AtomCoord:
    """3D atomic coordinate."""
    x: float
    y: float
    z: float

    def distance_to(self, other: 'AtomCoord') -> float:
        """Euclidean distance to another atom."""
        return math.sqrt(
            (self.x - other.x) ** 2 +
            (self.y - other.y) ** 2 +
            (self.z - other.z) ** 2
        )


@dataclass
class Residue:
    """Single residue with coordinates and confidence."""
    chain_id: str
    res_id: int
    res_name: str
    atoms: Dict[str, AtomCoord] = field(default_factory=dict)
    plddt: Optional[float] = None

    @property
    def ca_coord(self) -> Optional[AtomCoord]:
        """Get CA atom coordinate (protein) or C1' (RNA)."""
        if 'CA' in self.atoms:
            return self.atoms['CA']
        if "C1'" in self.atoms:
            return self.atoms["C1'"]
        # Fallback to any atom
        if self.atoms:
            return next(iter(self.atoms.values()))
        return None


@dataclass
class AF3Confidences:
    """Parsed AF3 confidence metrics."""
    plddt: List[float]  # Per-residue pLDDT
    pae: List[List[float]]  # Predicted Aligned Error matrix
    chain_pair_pae: Dict[str, float]  # Min PAE between chain pairs
    ptm: float  # Predicted TM-score
    iptm: float  # Interface pTM
    ranking_score: float


@dataclass
class AF3Structure:
    """Parsed AF3 structure with residues and chains."""
    residues: List[Residue]
    chains: Dict[str, List[Residue]]
    confidences: Optional[AF3Confidences] = None

    def get_chain(self, chain_id: str) -> List[Residue]:
        """Get all residues for a chain."""
        return self.chains.get(chain_id, [])

    def get_interface_contacts(
        self,
        chain1: str,
        chain2: str,
        distance_threshold: float = 8.0
    ) -> List[Tuple[Residue, Residue, float]]:
        """
        Find contacts between two chains.

        Args:
            chain1: First chain ID
            chain2: Second chain ID
            distance_threshold: Max distance for contact (Å)

        Returns:
            List of (residue1, residue2, distance) tuples
        """
        contacts = []
        residues1 = self.get_chain(chain1)
        residues2 = self.get_chain(chain2)

        for r1 in residues1:
            coord1 = r1.ca_coord
            if not coord1:
                continue

            for r2 in residues2:
                coord2 = r2.ca_coord
                if not coord2:
                    continue

                dist = coord1.distance_to(coord2)
                if dist <= distance_threshold:
                    contacts.append((r1, r2, dist))

        return contacts


class AF3Parser:
    """Parser for AlphaFold3 output files."""

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)

    def parse(self) -> Optional[AF3Structure]:
        """
        Parse AF3 outputs from directory.

        Returns:
            AF3Structure with residues and confidences
        """
        # Find structure file
        cif_files = list(self.output_dir.glob('*_model.cif'))
        if not cif_files:
            cif_files = list(self.output_dir.glob('*.cif'))

        if not cif_files:
            print(f"No CIF file found in {self.output_dir}", file=sys.stderr)
            return None

        structure = self._parse_cif(cif_files[0])

        # Parse confidences
        conf_files = list(self.output_dir.glob('*_confidences.json'))
        if not conf_files:
            conf_files = list(self.output_dir.glob('*confidence*.json'))

        if conf_files:
            structure.confidences = self._parse_confidences(conf_files[0])

            # Assign pLDDT to residues
            if structure.confidences and structure.confidences.plddt:
                for i, res in enumerate(structure.residues):
                    if i < len(structure.confidences.plddt):
                        res.plddt = structure.confidences.plddt[i]

        return structure

    def _parse_cif(self, cif_path: Path) -> AF3Structure:
        """Parse mmCIF file for coordinates."""
        residues = []
        chains: Dict[str, List[Residue]] = {}
        current_residue: Optional[Residue] = None

        with open(cif_path, 'r') as f:
            in_atom_site = False
            headers = []

            for line in f:
                line = line.strip()

                # Detect atom_site loop
                if line.startswith('loop_'):
                    in_atom_site = False
                    headers = []
                elif line.startswith('_atom_site.'):
                    in_atom_site = True
                    headers.append(line.split('.')[1])
                elif in_atom_site and line and not line.startswith('_') and not line.startswith('#'):
                    # Parse atom record
                    values = line.split()
                    if len(values) < len(headers):
                        continue

                    record = dict(zip(headers, values))

                    # Skip non-ATOM records
                    group = record.get('group_PDB', 'ATOM')
                    if group not in ['ATOM', 'HETATM']:
                        continue

                    chain_id = record.get('auth_asym_id', record.get('label_asym_id', 'A'))
                    res_id = int(record.get('auth_seq_id', record.get('label_seq_id', 0)))
                    res_name = record.get('auth_comp_id', record.get('label_comp_id', 'UNK'))
                    atom_name = record.get('auth_atom_id', record.get('label_atom_id', 'CA'))

                    try:
                        x = float(record.get('Cartn_x', 0))
                        y = float(record.get('Cartn_y', 0))
                        z = float(record.get('Cartn_z', 0))
                    except ValueError:
                        continue

                    # Create or update residue
                    if (current_residue is None or
                        current_residue.chain_id != chain_id or
                        current_residue.res_id != res_id):

                        if current_residue:
                            residues.append(current_residue)
                            if current_residue.chain_id not in chains:
                                chains[current_residue.chain_id] = []
                            chains[current_residue.chain_id].append(current_residue)

                        current_residue = Residue(
                            chain_id=chain_id,
                            res_id=res_id,
                            res_name=res_name
                        )

                    current_residue.atoms[atom_name] = AtomCoord(x, y, z)

            # Add last residue
            if current_residue:
                residues.append(current_residue)
                if current_residue.chain_id not in chains:
                    chains[current_residue.chain_id] = []
                chains[current_residue.chain_id].append(current_residue)

        return AF3Structure(residues=residues, chains=chains)

    def _parse_confidences(self, conf_path: Path) -> AF3Confidences:
        """Parse AF3 confidence JSON file."""
        with open(conf_path, 'r') as f:
            data = json.load(f)

        # Extract metrics
        plddt = data.get('atom_plddts', data.get('plddt', []))

        # PAE matrix
        pae = data.get('pae', [])

        # Chain pair PAE (min PAE between chains)
        chain_pair_pae = {}
        if 'chain_pair_pae_min' in data:
            chain_pair_pae = data['chain_pair_pae_min']
        elif pae:
            # Calculate from PAE matrix if chain info available
            pass

        ptm = data.get('ptm', 0.0)
        iptm = data.get('iptm', 0.0)
        ranking = data.get('ranking_score', data.get('ranking_confidence', 0.0))

        return AF3Confidences(
            plddt=plddt,
            pae=pae,
            chain_pair_pae=chain_pair_pae,
            ptm=ptm,
            iptm=iptm,
            ranking_score=ranking
        )


@dataclass
class BindingAnalysis:
    """Analysis of RNA-protein binding from AF3 prediction."""
    rna_chain: str
    protein_chain: str
    n_contacts: int
    min_contact_distance: float
    mean_contact_distance: float
    interface_plddt_rna: float
    interface_plddt_protein: float
    chain_pair_pae_min: float
    contact_residues_rna: List[int]
    contact_residues_protein: List[int]


def analyze_binding(
    structure: AF3Structure,
    rna_chain: str = "R",
    protein_chain: str = "P",
    contact_threshold: float = 8.0
) -> Optional[BindingAnalysis]:
    """
    Analyze RNA-protein binding from AF3 structure.

    Args:
        structure: Parsed AF3 structure
        rna_chain: RNA chain ID
        protein_chain: Protein chain ID
        contact_threshold: Distance threshold for contacts (Å)

    Returns:
        BindingAnalysis with binding metrics
    """
    # Get contacts
    contacts = structure.get_interface_contacts(
        rna_chain, protein_chain, contact_threshold
    )

    if not contacts:
        return None

    # Calculate metrics
    distances = [c[2] for c in contacts]
    rna_residues = list(set(c[0].res_id for c in contacts))
    protein_residues = list(set(c[1].res_id for c in contacts))

    # Interface pLDDT
    rna_chain_residues = structure.get_chain(rna_chain)
    protein_chain_residues = structure.get_chain(protein_chain)

    interface_rna = [r for r in rna_chain_residues if r.res_id in rna_residues]
    interface_protein = [r for r in protein_chain_residues if r.res_id in protein_residues]

    rna_plddt = [r.plddt for r in interface_rna if r.plddt is not None]
    protein_plddt = [r.plddt for r in interface_protein if r.plddt is not None]

    # Chain pair PAE
    chain_pair_pae = 0.0
    if structure.confidences and structure.confidences.chain_pair_pae:
        key = f"{rna_chain}_{protein_chain}"
        alt_key = f"{protein_chain}_{rna_chain}"
        chain_pair_pae = structure.confidences.chain_pair_pae.get(
            key, structure.confidences.chain_pair_pae.get(alt_key, 0.0)
        )

    return BindingAnalysis(
        rna_chain=rna_chain,
        protein_chain=protein_chain,
        n_contacts=len(contacts),
        min_contact_distance=min(distances),
        mean_contact_distance=sum(distances) / len(distances),
        interface_plddt_rna=sum(rna_plddt) / len(rna_plddt) if rna_plddt else 0.0,
        interface_plddt_protein=sum(protein_plddt) / len(protein_plddt) if protein_plddt else 0.0,
        chain_pair_pae_min=chain_pair_pae,
        contact_residues_rna=sorted(rna_residues),
        contact_residues_protein=sorted(protein_residues)
    )


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Parse AF3 outputs')
    parser.add_argument('output_dir', help='AF3 output directory')
    parser.add_argument('--rna-chain', default='R', help='RNA chain ID')
    parser.add_argument('--protein-chain', default='P', help='Protein chain ID')

    args = parser.parse_args()

    af3_parser = AF3Parser(args.output_dir)
    structure = af3_parser.parse()

    if structure:
        print(f"Parsed {len(structure.residues)} residues")
        print(f"Chains: {list(structure.chains.keys())}")

        if structure.confidences:
            print(f"pTM: {structure.confidences.ptm:.3f}")
            print(f"ipTM: {structure.confidences.iptm:.3f}")

        binding = analyze_binding(
            structure,
            rna_chain=args.rna_chain,
            protein_chain=args.protein_chain
        )

        if binding:
            print(f"\nBinding Analysis:")
            print(f"  Contacts: {binding.n_contacts}")
            print(f"  Min distance: {binding.min_contact_distance:.2f} Å")
            print(f"  Interface pLDDT (RNA): {binding.interface_plddt_rna:.1f}")
            print(f"  Interface pLDDT (protein): {binding.interface_plddt_protein:.1f}")
            print(f"  Chain pair PAE min: {binding.chain_pair_pae_min:.2f}")
