import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def feature_hydrophobicity(protein):
    """Extract hydrophobicity-based features"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Create a feature vector with hydrophobicity-related metrics
        features = [
            analysis.gravy(),  # Grand average of hydropathy
            analysis.aromaticity(),
            analysis.instability_index(),
            analysis.isoelectric_point(),
            analysis.charge_at_pH(7.0),
        ]
        # Pad to 10 features
        features.extend([0.0] * (10 - len(features)))
        return np.array(features[:10])
    except:
        return np.zeros(10)


def feature_molecular_weight(protein):
    """Extract molecular weight and composition features"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Get amino acid percentages for common residues
        aa_percent = analysis.get_amino_acids_percent()
        features = [
            analysis.molecular_weight() / 10000.0,  # Normalize
            aa_percent.get('A', 0.0),  # Alanine
            aa_percent.get('L', 0.0),  # Leucine
            aa_percent.get('G', 0.0),  # Glycine
            aa_percent.get('V', 0.0),  # Valine
            aa_percent.get('E', 0.0),  # Glutamate
            aa_percent.get('K', 0.0),  # Lysine
            aa_percent.get('D', 0.0),  # Aspartate
            aa_percent.get('S', 0.0),  # Serine
            aa_percent.get('T', 0.0),  # Threonine
        ]
        return np.array(features[:10])
    except:
        return np.zeros(10)


def feature_secondary_structure(protein):
    """Estimate secondary structure propensities"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Get secondary structure fractions
        helix, turn, sheet = analysis.secondary_structure_fraction()
        
        # Additional structural features
        features = [
            helix,
            turn,
            sheet,
            len(seq) / 100.0,  # Normalized length
            analysis.molar_extinction_coefficient()[0] / 10000.0,
            analysis.molar_extinction_coefficient()[1] / 10000.0,
            0.0, 0.0, 0.0, 0.0  # Padding
        ]
        return np.array(features[:10])
    except:
        return np.zeros(10)


FEATURE_FUNCTIONS = [
    feature_hydrophobicity,
    feature_molecular_weight,
    feature_secondary_structure
]
