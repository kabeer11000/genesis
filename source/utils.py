from Bio.PDB import PDBParser, PPBuilder
import numpy as np


def load_proteins(paths):
    parser = PDBParser(QUIET=True)
    proteins = []
    for path in paths:
        structure = parser.get_structure('protein', path)
        model = structure[0]
        chain = next(model.get_chains())
        
        # Get CA atoms and coordinates
        ca_atoms = []
        coords = []
        for residue in chain:
            if 'CA' in residue:
                ca_atoms.append(residue['CA'])
                coords.append(residue['CA'].get_coord())
        
        # Get sequence using PPBuilder for single-letter codes
        ppb = PPBuilder()
        peptides = ppb.build_peptides(chain)
        seq = ""
        for peptide in peptides:
            seq += str(peptide.get_sequence())
        
        proteins.append({
            'structure': structure,
            'ca_atoms': ca_atoms,
            'coords': np.array(coords),
            'sequence': seq
        })
    return proteins


def rmsd(pred_coords, true_coords):
    """Calculate RMSD between two sets of coordinates"""
    if pred_coords.shape != true_coords.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    # Center both structures
    pred_centered = pred_coords - np.mean(pred_coords, axis=0)
    true_centered = true_coords - np.mean(true_coords, axis=0)
    
    # Calculate RMSD
    diff = pred_centered - true_centered
    rmsd_value = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd_value
