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
        
        # Enhanced common drugs database as fallback for cloud deployment
        common_drugs = {
            # Pain relievers
            "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 
            "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
            "paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
            "naproxen": "COC1=CC=CC(=C1)C(C)C(=O)O",
            
            # Blood thinners
            "warfarin": "CC1=C(C2=C(C=CC=C2OC1=O)O)C(CC(C3=CC=CC=C3)C4=CC=CC=C4)=O",
            "heparin": "CC(=O)NC1C(C(C(OC1O)COS(=O)(=O)O)O)O",
            
            # Diabetes
            "metformin": "CN(C)C(=N)NC(=N)N",
            "insulin": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
            
            # Stimulants
            "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            
            # Simple molecules
            "ethanol": "CCO",
            "glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
            "water": "O",
            "methanol": "CO",
            
            # Antibiotics
            "penicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
            "amoxicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
            
            # Cardiovascular
            "atenolol": "CC(C)NCC(COC1=CC=C(C=C1)CC(=O)N)O",
            "lisinopril": "CCCCN1CCCC1C(=O)N2CCCC2C(=O)O",
            
            # Mental health
            "fluoxetine": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
            "sertraline": "CNC1CCC(C2=CC=C(Cl)C=C12)C3=CC=C(Cl)C=C3"
        }
        
        # Try to find the drug in our fallback database
        drug_name_lower = smiles_or_name.lower().strip()
        if drug_name_lower in common_drugs:
            try:
                smiles = common_drugs[drug_name_lower]
                print(f"Using fallback SMILES for {smiles_or_name}: {smiles}")
                
                if RDKIT_AVAILABLE:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        return {
                            "smiles": smiles,
                            "name": smiles_or_name.title(),
                            "mw": Descriptors.MolWt(mol),
                            "logp": Descriptors.MolLogP(mol),
                            "hbd": Descriptors.NumHDonors(mol),
                            "hba": Descriptors.NumHAcceptors(mol),
                            "image": Draw.MolToImage(mol, size=(300, 300))
                        }
                    else:
                        return _estimate_molecular_properties_fallback(smiles)
                else:
                    # Use fallback estimation
                    result = _estimate_molecular_properties_fallback(smiles)
                    if result:
                        result['name'] = smiles_or_name.title()
                    return result
            except Exception as fallback_error:
                print(f"Fallback processing failed for {smiles_or_name}: {fallback_error}")
        
        return None