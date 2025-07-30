import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from molecule_utils import get_molecule_features
from quantum_encoder import encode_features
from interaction_score import compute_fidelity

def interpret_score(fidelity):
    if fidelity > 0.85:
        return "🔴 High Risk of Interaction"
    elif fidelity > 0.6:
        return "🟠 Moderate Risk"
    else:
        return "🟢 Low or No Risk"

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol) if mol else None

st.set_page_config(page_title="Q-Cure | Quantum Drug Interaction Predictor")

st.title("🧪 Q-Cure: Quantum-Powered Drug Interaction Predictor")
st.markdown("Predict potential interactions between two drugs using quantum state similarity.")

# Input section
col1, col2 = st.columns(2)
with col1:
    smiles1 = st.text_input("Enter SMILES for Drug 1 (e.g. CCO):", "CCO")
with col2:
    smiles2 = st.text_input("Enter SMILES for Drug 2 (e.g. CC(=O)O):", "CC(=O)O")

if st.button("🔍 Predict Interaction"):
    features1 = get_molecule_features(smiles1)
    features2 = get_molecule_features(smiles2)

    if features1 and features2:
        _, sv1 = encode_features(features1)
        _, sv2 = encode_features(features2)
        fidelity = compute_fidelity(sv1, sv2)
        interpretation = interpret_score(fidelity)

        st.success(f"**Quantum Fidelity Score**: `{fidelity:.4f}`")
        st.info(f"**Risk Assessment**: {interpretation}")

        col1, col2 = st.columns(2)
        with col1:
            st.image(draw_molecule(smiles1), caption="Drug 1 Structure")
        with col2:
            st.image(draw_molecule(smiles2), caption="Drug 2 Structure")

    else:
        st.error("❌ Invalid SMILES input. Please check your entries.")
