from rdkit import Chem
from rdkit.Chem import Descriptors

def get_molecule_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mw = Descriptors.MolWt(mol)             # Molecular weight
    logp = Descriptors.MolLogP(mol)         # Hydrophobicity
    hbd = Descriptors.NumHDonors(mol)       # Hydrogen bond donors
    hba = Descriptors.NumHAcceptors(mol)    # Hydrogen bond acceptors

    return [mw, logp, hbd, hba]
