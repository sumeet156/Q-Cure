from molecule_utils import get_molecule_features
from quantum_encoder import encode_features
from interaction_score import compute_fidelity

def interpret_score(fidelity):
    if fidelity > 0.85:
        return "High Risk of Interaction"
    elif fidelity > 0.6:
        return "Moderate Risk"
    else:
        return "Low or No Risk"

if __name__ == "__main__":
    # Test multiple molecule pairs
    test_pairs = [
        ("CCO", "CCO"),           # Same molecule (ethanol vs ethanol)
        ("CCO", "CC(=O)O"),       # Different molecules (ethanol vs acetic acid)
        ("CCO", "CCC"),           # Ethanol vs propane
        ("CC(=O)O", "CC(=O)O"),   # Same molecule (acetic acid vs acetic acid)
    ]
    
    for smiles1, smiles2 in test_pairs:
        features1 = get_molecule_features(smiles1)
        features2 = get_molecule_features(smiles2)

        if features1 and features2:
            _, sv1 = encode_features(features1)
            _, sv2 = encode_features(features2)

            fidelity = compute_fidelity(sv1, sv2)
            interpretation = interpret_score(fidelity)
            print(f"{smiles1} vs {smiles2}: Fidelity = {fidelity:.4f} ({interpretation})")
        else:
            print(f"Invalid SMILES: {smiles1} or {smiles2}")
