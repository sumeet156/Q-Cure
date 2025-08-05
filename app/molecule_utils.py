# Handle RDKit import gracefully for cloud deployment
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw
    RDKIT_AVAILABLE = True
except ImportError:
    print("RDKit not available - using simplified molecular processing")
    RDKIT_AVAILABLE = False
    Chem = None
    Descriptors = None
    Draw = None

import pubchempy as pcp
import time
import io
import base64

def _estimate_molecular_properties_fallback(smiles):
    """Estimate molecular properties when RDKit is not available."""
    if not smiles:
        return None
    
    # Simple heuristics based on SMILES string
    carbon_count = smiles.count('C') + smiles.count('c')
    nitrogen_count = smiles.count('N') + smiles.count('n')
    oxygen_count = smiles.count('O') + smiles.count('o')
    
    # Estimate molecular weight
    estimated_mw = (carbon_count * 12) + (nitrogen_count * 14) + (oxygen_count * 16)
    estimated_mw += len(smiles) * 1  # Add for other atoms
    estimated_mw = max(estimated_mw, 50)
    
    # Estimate other properties
    estimated_logp = (carbon_count - oxygen_count * 2) / max(carbon_count + oxygen_count, 1)
    estimated_logp = max(-5, min(5, estimated_logp))
    
    estimated_hbd = smiles.count('OH') + smiles.count('NH')
    estimated_hba = oxygen_count + nitrogen_count
    
    return {
        'smiles': smiles,
        'mw': estimated_mw,
        'logp': estimated_logp,
        'hbd': estimated_hbd,
        'hba': estimated_hba,
        'name': 'Compound (estimated)',
        'image': None
    }

def get_molecule_features(smiles_or_name: str):
    """Fetch molecule properties from SMILES or drug name."""
    if not smiles_or_name or not smiles_or_name.strip():
        return None
        
    smiles_or_name = smiles_or_name.strip()
    
    try:
        # Enhanced SMILES detection
        smiles_chars = ["=", "(", ")", "[", "]", "#", "@", "+", "-"]
        is_likely_smiles = any(c in smiles_or_name for c in smiles_chars) or smiles_or_name.isupper()
        
        if not is_likely_smiles:
            # Try to fetch from PubChem by name
            print(f"Searching PubChem for: {smiles_or_name}")
            try:
                compounds = pcp.get_compounds(smiles_or_name, "name")
                
                if not compounds:
                    # Try alternative search methods
                    compounds = pcp.get_compounds(smiles_or_name, "name", searchtype="similarity")
                    
                if compounds and len(compounds) > 0:
                    smiles = compounds[0].canonical_smiles
                    if smiles:
                        print(f"Found SMILES: {smiles}")
                    else:
                        print(f"Found compound but no SMILES available")
                        raise Exception("No SMILES found")
                else:
                    print(f"No compounds found for: {smiles_or_name}")
                    raise Exception("No compounds found")
                    
            except Exception as e:
                print(f"PubChem search failed: {e}")
                # Fall back to common drugs dictionary
                raise Exception("PubChem lookup failed")
        else:
            smiles = smiles_or_name
        
        # Process molecule based on available libraries
        if RDKIT_AVAILABLE:
            # Validate and create molecule using RDKit
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                properties = {
                    "smiles": smiles,
                    "mw": Descriptors.MolWt(mol),
                    "logp": Descriptors.MolLogP(mol),
                    "hbd": Descriptors.NumHDonors(mol),
                    "hba": Descriptors.NumHAcceptors(mol),
                    "image": Draw.MolToImage(mol, size=(300, 300))
                }
            else:
                print(f"RDKit failed to parse SMILES: {smiles}")
                # Fall back to estimation
                properties = _estimate_molecular_properties_fallback(smiles)
        else:
            # Use fallback estimation when RDKit is not available
            properties = _estimate_molecular_properties_fallback(smiles)
        
        if properties:
            print(f"Successfully processed: {smiles_or_name} -> MW: {properties['mw']:.2f}")
        return properties
        
    except Exception as e:
        print(f"Error fetching molecule '{smiles_or_name}': {e}")
        
        # Try basic common drugs as fallback
        common_drugs = {
            "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 
            "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
            "paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
            "warfarin": "CC1=C(C2=C(C=CC=C2OC1=O)O)C(=O)CC(C3=CC=CC=C3)C4=CC=CC=C4",
            "metformin": "CN(C)C(=N)N=C(N)N",
            "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "ethanol": "CCO",
            "glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
        }
        
        if smiles_or_name.lower() in common_drugs:
            try:
                smiles = common_drugs[smiles_or_name.lower()]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return {
                        "smiles": smiles,
                        "mw": Descriptors.MolWt(mol),
                        "logp": Descriptors.MolLogP(mol),
                        "hbd": Descriptors.NumHDonors(mol),
                        "hba": Descriptors.NumHAcceptors(mol),
                        "image": Draw.MolToImage(mol, size=(300, 300))
                    }
            except:
                pass
        
        return None