# 🧪 Q-Cure: Quantum-Powered Drug Interaction Predictor

Q-Cure is a hybrid quantum-classical platform that predicts potential **adverse drug interactions** (ADIs) using molecular features encoded into quantum states. The project leverages **Qiskit**, **RDKit**, and **Streamlit** to deliver an interactive, educational, and extensible proof-of-concept in the field of **quantum-enhanced chemoinformatics**.

---

## 🚀 Overview

Adverse drug interactions are responsible for thousands of hospitalizations and deaths annually. Existing computational tools struggle with scalability and predictive precision. **Q-Cure** explores how **quantum computing** can assist in detecting molecular conflicts by simulating **quantum similarity** between drug representations.

---

## 🧠 How It Works

### 🔬 Molecular Processing (RDKit)
- Users input **SMILES strings** for two drugs.
- RDKit extracts key descriptors:
  - Molecular Weight
  - LogP (hydrophobicity)
  - Hydrogen Bond Donors / Acceptors

### 🔗 Quantum Encoding (Qiskit)
- Descriptors are embedded into a **quantum circuit** using a `ZZFeatureMap`.
- Circuits are simulated on a **statevector simulator**.
- The **fidelity** between two quantum statevectors is used to assess interaction risk.

### 🧠 Interpretation Logic
- **Fidelity > 0.85** → 🔴 High Risk
- **Fidelity > 0.6** → 🟠 Moderate Risk
- **Fidelity ≤ 0.6** → 🟢 Low/No Risk

---


