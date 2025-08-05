import streamlit as st
import sys
import os

# Add the app directory to the Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from quantum_chemistry import QuantumChemistry
from molecule_utils import get_molecule_features
from interaction_score import assess_risk

# Page configuration
st.set_page_config(
    page_title="Q-Cure | Quantum Drug Interaction Predictor",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better UI
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
    }
    .stButton > button {
        background-color: #2E86AB;
        color: white;
        border: none;
        border-radius: 10px;
        padding: 0.5rem 1rem;
        font-weight: bold;
    }
    .stButton > button:hover {
        background-color: #1F5F8B;
    }
    .feature-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 10px;
        border-left: 4px solid #2E86AB;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Main title
st.markdown('<h1 class="main-header">🧪 Q-Cure: Quantum-Powered Drug Interaction Predictor</h1>', unsafe_allow_html=True)

# Sidebar information
with st.sidebar:
    st.header("⚙️ Configuration")
    
    # IBM Quantum toggle
    use_ibmq = st.checkbox("🚀 Use IBM Quantum Hardware", False, 
                          help="Enable IBM Quantum hardware for real quantum computation (slower but more accurate)")
    
    if use_ibmq:
        st.info("🔬 IBM Quantum mode enabled. Computations will use real quantum hardware when available.")
    else:
        st.info("💻 Using quantum simulator for fast computations.")
    
    st.header("📚 About Q-Cure")
    st.markdown("""
    **Q-Cure** leverages quantum computing to predict drug-drug interactions:
    
    ✅ **Quantum Chemistry (VQE)** – Simulates molecular energies using quantum circuits
    
    ✅ **Real Drug Data (PubChem)** – Fetches real molecular structures and properties
    
    ✅ **IBM Quantum Hardware** – Optional real quantum execution
    
    ✅ **Risk Assessment** – Clear visual output with color-coded risk levels
    """)

# Initialize quantum backend
@st.cache_resource
def init_quantum_chemistry(use_ibmq):
    return QuantumChemistry(use_ibmq=use_ibmq)

qc = init_quantum_chemistry(use_ibmq)

# Main interface
st.header("🔬 Drug Interaction Analysis")

# Input section
col1, col2 = st.columns(2)

# Get values from session state if available, otherwise use defaults
default_drug1 = st.session_state.get("example_drug1", "aspirin")
default_drug2 = st.session_state.get("example_drug2", "ibuprofen")

with col1:
    st.subheader("💊 Drug 1")
    drug1 = st.text_input(
        "Enter drug name or SMILES:",
        value=default_drug1,
        placeholder="e.g., aspirin, ibuprofen, CCO",
        key="drug1_input"
    )
    
with col2:
    st.subheader("💊 Drug 2") 
    drug2 = st.text_input(
        "Enter drug name or SMILES:",
        value=default_drug2, 
        placeholder="e.g., warfarin, C1=CC=CC=C1",
        key="drug2_input"
    )

# Example drugs
st.markdown("**💡 Try these examples:**")
example_col1, example_col2, example_col3 = st.columns(3)

with example_col1:
    if st.button("Aspirin + Warfarin", key="example1"):
        st.session_state["example_drug1"] = "aspirin"
        st.session_state["example_drug2"] = "warfarin"
        st.rerun()

with example_col2:
    if st.button("Ibuprofen + Acetaminophen", key="example2"):
        st.session_state["example_drug1"] = "ibuprofen"
        st.session_state["example_drug2"] = "acetaminophen"
        st.rerun()

with example_col3:
    if st.button("Metformin + Caffeine", key="example3"):
        st.session_state["example_drug1"] = "metformin"
        st.session_state["example_drug2"] = "caffeine"
        st.rerun()

# Analysis button
if st.button("🔍 Predict Interaction", type="primary"):
    if drug1 and drug2:
        # Fetch molecular data
        with st.spinner("📊 Fetching molecular data from PubChem..."):
            mol1 = get_molecule_features(drug1)
            mol2 = get_molecule_features(drug2)
        
        if mol1 and mol2:
            # Display molecular info
            st.subheader("🧬 Molecular Information")
            
            mol_col1, mol_col2 = st.columns(2)
            
            with mol_col1:
                st.markdown(f"**{drug1.title()}**")
                st.markdown(f"- **SMILES:** `{mol1['smiles']}`")
                st.markdown(f"- **Molecular Weight:** {mol1['mw']:.2f} g/mol")
                st.markdown(f"- **LogP:** {mol1['logp']:.2f}")
                st.markdown(f"- **H-Bond Donors:** {mol1['hbd']}")
                st.markdown(f"- **H-Bond Acceptors:** {mol1['hba']}")
                
                # Display molecular structure image only if available
                if mol1.get("image") is not None:
                    st.image(mol1["image"], caption=f"{drug1.title()} Structure", width=300)
                else:
                    st.info(f"📋 {drug1.title()} molecular structure\n\nSMILES: `{mol1['smiles']}`")
            
            with mol_col2:
                st.markdown(f"**{drug2.title()}**")
                st.markdown(f"- **SMILES:** `{mol2['smiles']}`")
                st.markdown(f"- **Molecular Weight:** {mol2['mw']:.2f} g/mol")
                st.markdown(f"- **LogP:** {mol2['logp']:.2f}")
                st.markdown(f"- **H-Bond Donors:** {mol2['hbd']}")
                st.markdown(f"- **H-Bond Acceptors:** {mol2['hba']}")
                
                # Display molecular structure image only if available
                if mol2.get("image") is not None:
                    st.image(mol2["image"], caption=f"{drug2.title()} Structure", width=300)
                else:
                    st.info(f"📋 {drug2.title()} molecular structure\n\nSMILES: `{mol2['smiles']}`")
            
            # Quantum computation
            st.subheader("⚛️ Quantum Chemistry Analysis")
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            with st.spinner("🔬 Running quantum molecular simulations..."):
                status_text.text("Computing molecular energies...")
                progress_bar.progress(25)
                
                result = qc.compute_interaction(mol1["smiles"], mol2["smiles"])
                progress_bar.progress(75)
                
                status_text.text("Assessing interaction risk...")
                risk = assess_risk(result["interaction_energy"])
                progress_bar.progress(100)
                
                status_text.text("Analysis complete!")
            
            # Results
            st.subheader("📊 Results")
            
            # Energy results
            energy_col1, energy_col2, energy_col3 = st.columns(3)
            
            with energy_col1:
                st.metric(
                    label="Drug 1 Energy",
                    value=f"{result['energy1']:.4f} Ha",
                    help="Molecular energy of the first drug"
                )
            
            with energy_col2:
                st.metric(
                    label="Drug 2 Energy", 
                    value=f"{result['energy2']:.4f} Ha",
                    help="Molecular energy of the second drug"
                )
            
            with energy_col3:
                st.metric(
                    label="Interaction Energy",
                    value=f"{result['interaction_energy']:.4f} Ha",
                    delta=f"Complex: {result['energy_complex']:.4f} Ha",
                    help="Energy change upon drug interaction"
                )
            
            # Risk assessment
            st.subheader("⚠️ Risk Assessment")
            
            if "High Risk" in risk:
                st.error(f"**{risk}**")
                st.markdown("""
                **⚠️ High Risk Interaction Detected!**
                - Strong binding energy indicates potential drug interaction
                - Consult healthcare provider before combining these medications
                - Monitor for side effects if used together
                """)
            elif "Moderate Risk" in risk:
                st.warning(f"**{risk}**")
                st.markdown("""
                **🟠 Moderate Risk Interaction**
                - Some interaction potential exists
                - Exercise caution when combining
                - Consider timing of administration
                """)
            else:
                st.success(f"**{risk}**")
                st.markdown("""
                **✅ Low Risk Interaction**
                - Minimal interaction expected
                - Generally safe to combine
                - Follow standard dosing guidelines
                """)
            
            # Technical details
            with st.expander("🔬 Technical Details"):
                st.markdown(f"""
                **Quantum Computation Details:**
                - **Method:** Variational Quantum Eigensolver (VQE) simulation
                - **Backend:** {'IBM Quantum Hardware' if use_ibmq else 'Quantum Simulator'}
                - **Molecular Encoding:** 4-qubit quantum circuit per molecule
                - **Features Encoded:** Molecular weight, LogP, H-bond donors/acceptors
                
                **Energy Components (Hartree units):**
                - Individual Drug Energies: {result['energy1']:.6f}, {result['energy2']:.6f}
                - Complex Formation Energy: {result['energy_complex']:.6f}
                - Net Interaction Energy: {result['interaction_energy']:.6f}
                
                **Interpretation:**
                - Negative interaction energy indicates favorable binding
                - More negative values suggest stronger interactions
                - Values below -0.5 Ha indicate high interaction risk
                """)
        
        else:
            st.error("❌ Could not fetch molecular data. Please check drug names or SMILES notation.")
            st.info("💡 **Tip:** Try common drug names like 'aspirin', 'ibuprofen', or valid SMILES strings like 'CCO' for ethanol.")
    
    else:
        st.warning("⚠️ Please enter both drugs to analyze interactions.")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>🧪 <strong>Q-Cure</strong> | Quantum-Powered Drug Interaction Predictor</p>
    <p><em>Disclaimer: This tool is for research purposes only. Always consult healthcare professionals for medical advice.</em></p>
</div>
""", unsafe_allow_html=True)
