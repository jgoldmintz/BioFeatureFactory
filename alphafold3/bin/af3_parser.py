# BioFeatureFactory
# Copyright (C) 2023-2026  Jacob Goldmintz
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
import statistics
from collections import Counter


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
            distance_threshold: Max distance for contact (A)

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
        # Find structure file (AF3 writes into a timestamped subdirectory)
        cif_files = list(self.output_dir.glob('*_model.cif'))
        if not cif_files:
            cif_files = list(self.output_dir.glob('**/*_model.cif'))
        # Filter out per-sample files, prefer the top-level ranked model
        top_level = [f for f in cif_files if 'seed-' not in f.name]
        if top_level:
            cif_files = top_level

        if not cif_files:
            print(f"No CIF file found in {self.output_dir}", file=sys.stderr)
            return None

        structure = self._parse_cif(cif_files[0])

        # Parse confidences from both files (same directory as the CIF)
        cif_parent = cif_files[0].parent

        # Per-atom data (_confidences.json): plddt, pae matrix
        conf_files = list(cif_parent.glob('*_confidences.json'))
        if not conf_files:
            conf_files = list(self.output_dir.glob('**/*_confidences.json'))
            conf_files = [f for f in conf_files if 'seed-' not in f.name and 'summary' not in f.name] or conf_files
        # Exclude summary files from per-atom list
        conf_files = [f for f in conf_files if 'summary' not in f.name]

        # Summary data (_summary_confidences.json): ptm, iptm, chain_pair_pae_min
        summary_files = list(cif_parent.glob('*_summary_confidences.json'))
        if not summary_files:
            summary_files = list(self.output_dir.glob('**/*_summary_confidences.json'))
            summary_files = [f for f in summary_files if 'seed-' not in f.name] or summary_files

        if conf_files or summary_files:
            structure.confidences = self._parse_confidences(
                conf_files[0] if conf_files else None,
                summary_files[0] if summary_files else None
            )

            # Assign pLDDT to residues (atom_plddts is per-atom; average per residue)
            if structure.confidences and structure.confidences.plddt:
                atom_idx = 0
                for res in structure.residues:
                    n_atoms = len(res.atoms)
                    if atom_idx + n_atoms <= len(structure.confidences.plddt):
                        res_plddts = structure.confidences.plddt[atom_idx:atom_idx + n_atoms]
                        res.plddt = sum(res_plddts) / len(res_plddts)
                    atom_idx += n_atoms

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

    def _parse_confidences(self, conf_path: Optional[Path], summary_path: Optional[Path] = None) -> AF3Confidences:
        """Parse AF3 confidence files.

        Args:
            conf_path: Per-atom confidences (_confidences.json) — plddt, pae matrix
            summary_path: Summary confidences (_summary_confidences.json) — ptm, iptm, chain_pair_pae_min
        """
        plddt = []
        pae = []

        if conf_path:
            with open(conf_path, 'r') as f:
                data = json.load(f)
            plddt = data.get('atom_plddts', data.get('plddt', []))
            pae = data.get('pae', [])

        # Aggregate metrics from summary file
        chain_pair_pae = {}
        ptm = 0.0
        iptm = 0.0
        ranking = 0.0

        if summary_path:
            with open(summary_path, 'r') as f:
                summary = json.load(f)

            ptm = summary.get('ptm', 0.0)
            iptm = summary.get('iptm', 0.0)
            ranking = summary.get('ranking_score', 0.0)

            if 'chain_pair_pae_min' in summary:
                raw = summary['chain_pair_pae_min']
                if isinstance(raw, dict):
                    chain_pair_pae = raw
                elif isinstance(raw, list) and len(raw) >= 2:
                    # Matrix format: row=chain_i, col=chain_j
                    # For a 2-chain complex (RNA=0, Protein=1), cross-chain values
                    chain_pair_pae = {
                        'cross_chain_min': min(raw[0][1] if len(raw[0]) > 1 else 999,
                                              raw[1][0] if len(raw[1]) > 0 else 999)
                    }

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
        contact_threshold: Distance threshold for contacts (A)

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
        cpae = structure.confidences.chain_pair_pae
        key = f"{rna_chain}_{protein_chain}"
        alt_key = f"{protein_chain}_{rna_chain}"
        chain_pair_pae = cpae.get(key, cpae.get(alt_key, cpae.get('cross_chain_min', 0.0)))

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


@dataclass
class AggregatedBindingAnalysis:
    """Aggregated binding analysis across multiple AF3 samples."""
    mean: BindingAnalysis
    std_n_contacts: float
    std_min_contact_distance: float
    std_interface_plddt_rna: float
    std_interface_plddt_protein: float
    std_chain_pair_pae_min: float
    n_samples: int
    contact_frequency_rna: Dict[int, float]
    contact_frequency_protein: Dict[int, float]


def parse_all_samples(output_dir: str) -> List[AF3Structure]:
    """
    Parse all seed-N_sample-N subdirectories in an AF3 output directory.

    Falls back to top-ranked model if no sample subdirs exist.
    """
    base = Path(output_dir)
    sample_dirs = sorted(set(
        p.parent for p in base.glob('**/seed-*_sample-*/*_model.cif')
    ))

    if not sample_dirs:
        # Fall back to top-ranked model
        parser = AF3Parser(output_dir)
        structure = parser.parse()
        return [structure] if structure else []

    structures = []
    for sd in sample_dirs:
        parser = AF3Parser(str(sd))
        structure = parser.parse()
        if structure:
            structures.append(structure)

    return structures


def aggregate_binding_analyses(
    analyses: List[Optional[BindingAnalysis]],
    rna_chain: str = "R",
    protein_chain: str = "P"
) -> Optional[AggregatedBindingAnalysis]:
    """
    Aggregate multiple BindingAnalysis objects into mean/std.

    Filters out None entries (samples with no binding).
    Returns None if all samples have no binding.
    """
    valid = [a for a in analyses if a is not None]
    if not valid:
        return None

    n = len(valid)

    contacts_list = [a.n_contacts for a in valid]
    min_dist_list = [a.min_contact_distance for a in valid]
    plddt_rna_list = [a.interface_plddt_rna for a in valid]
    plddt_prot_list = [a.interface_plddt_protein for a in valid]
    pae_list = [a.chain_pair_pae_min for a in valid]

    # Contact frequency: fraction of samples each residue appears as contact
    rna_counter: Counter = Counter()
    prot_counter: Counter = Counter()
    for a in valid:
        rna_counter.update(a.contact_residues_rna)
        prot_counter.update(a.contact_residues_protein)

    contact_freq_rna = {rid: cnt / n for rid, cnt in rna_counter.items()}
    contact_freq_prot = {rid: cnt / n for rid, cnt in prot_counter.items()}

    # Union of all contact residues across samples
    all_rna = sorted(rna_counter.keys())
    all_prot = sorted(prot_counter.keys())

    mean_binding = BindingAnalysis(
        rna_chain=rna_chain,
        protein_chain=protein_chain,
        n_contacts=int(statistics.mean(contacts_list)),
        min_contact_distance=statistics.mean(min_dist_list),
        mean_contact_distance=statistics.mean([a.mean_contact_distance for a in valid]),
        interface_plddt_rna=statistics.mean(plddt_rna_list),
        interface_plddt_protein=statistics.mean(plddt_prot_list),
        chain_pair_pae_min=statistics.mean(pae_list),
        contact_residues_rna=all_rna,
        contact_residues_protein=all_prot
    )

    return AggregatedBindingAnalysis(
        mean=mean_binding,
        std_n_contacts=statistics.pstdev(contacts_list),
        std_min_contact_distance=statistics.pstdev(min_dist_list),
        std_interface_plddt_rna=statistics.pstdev(plddt_rna_list),
        std_interface_plddt_protein=statistics.pstdev(plddt_prot_list),
        std_chain_pair_pae_min=statistics.pstdev(pae_list),
        n_samples=n,
        contact_frequency_rna=contact_freq_rna,
        contact_frequency_protein=contact_freq_prot
    )


@dataclass
class InterfaceSite:
    """Per-residue data at or near the RNA-protein interface."""
    chain: str
    res_id: int
    res_name: str
    plddt: float
    is_contact: bool
    min_contact_distance: float


def extract_interface_sites(
    structure: AF3Structure,
    rna_chain: str = "R",
    protein_chain: str = "P",
    contact_threshold: float = 8.0,
    near_threshold: float = 12.0
) -> List[InterfaceSite]:
    """
    Extract per-residue data for residues at or near the interface.

    Includes all residues within near_threshold of the other chain.
    Marks is_contact=True for residues within contact_threshold.
    """
    rna_residues = structure.get_chain(rna_chain)
    prot_residues = structure.get_chain(protein_chain)
    sites = []

    # RNA residues: min distance to any protein residue
    for r in rna_residues:
        coord = r.ca_coord
        if not coord:
            continue
        min_dist = float('inf')
        for p in prot_residues:
            pc = p.ca_coord
            if pc:
                d = coord.distance_to(pc)
                if d < min_dist:
                    min_dist = d
        if min_dist <= near_threshold:
            sites.append(InterfaceSite(
                chain=rna_chain,
                res_id=r.res_id,
                res_name=r.res_name,
                plddt=r.plddt if r.plddt is not None else 0.0,
                is_contact=min_dist <= contact_threshold,
                min_contact_distance=round(min_dist, 2)
            ))

    # Protein residues: min distance to any RNA residue
    for p in prot_residues:
        coord = p.ca_coord
        if not coord:
            continue
        min_dist = float('inf')
        for r in rna_residues:
            rc = r.ca_coord
            if rc:
                d = coord.distance_to(rc)
                if d < min_dist:
                    min_dist = d
        if min_dist <= near_threshold:
            sites.append(InterfaceSite(
                chain=protein_chain,
                res_id=p.res_id,
                res_name=p.res_name,
                plddt=p.plddt if p.plddt is not None else 0.0,
                is_contact=min_dist <= contact_threshold,
                min_contact_distance=round(min_dist, 2)
            ))

    return sites


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
            print(f"  Min distance: {binding.min_contact_distance:.2f} A")
            print(f"  Interface pLDDT (RNA): {binding.interface_plddt_rna:.1f}")
            print(f"  Interface pLDDT (protein): {binding.interface_plddt_protein:.1f}")
            print(f"  Chain pair PAE min: {binding.chain_pair_pae_min:.2f}")
