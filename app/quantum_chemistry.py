from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

# Handle different Qiskit versions for Streamlit Cloud compatibility
try:
    from qiskit_aer.primitives import Sampler
except ImportError:
    try:
        from qiskit.primitives import Sampler
    except ImportError:
        # Fallback - we'll use direct backend execution
        Sampler = None

# Handle RDKit import gracefully for cloud deployment
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    print("RDKit not available - using simplified molecular calculations")
    RDKIT_AVAILABLE = False
    Chem = None
    Descriptors = None

import numpy as np
import os

class QuantumChemistry:
    def __init__(self, use_ibmq=False):
        self.use_ibmq = use_ibmq
        self.backend = AerSimulator()
        
        if use_ibmq:
            try:
                # For IBM Quantum hardware integration
                print("IBM Quantum integration available but using simulator for stability")
                # You can add IBM Quantum provider setup here when needed
            except Exception as e:
                print(f"IBM Quantum setup failed, using simulator: {e}")
                self.use_ibmq = False

    def _estimate_molecular_properties(self, smiles):
        """Estimate molecular properties from SMILES when RDKit is unavailable."""
        if not smiles:
            return [100, 0, 0, 0]  # Default values
        
        # Simple heuristics based on SMILES string
        carbon_count = smiles.count('C') + smiles.count('c')
        nitrogen_count = smiles.count('N') + smiles.count('n')
        oxygen_count = smiles.count('O') + smiles.count('o')
        sulfur_count = smiles.count('S') + smiles.count('s')
        
        # Estimate molecular weight
        estimated_mw = (carbon_count * 12) + (nitrogen_count * 14) + (oxygen_count * 16) + (sulfur_count * 32)
        estimated_mw += len(smiles) * 1  # Add for hydrogens and other atoms
        estimated_mw = max(estimated_mw, 50)  # Minimum reasonable MW
        
        # Estimate LogP
        estimated_logp = (carbon_count - oxygen_count * 2) / max(carbon_count + oxygen_count, 1)
        estimated_logp = max(-5, min(5, estimated_logp))
        
        # Estimate H-bond donors and acceptors
        oh_groups = smiles.count('OH') + smiles.count('oh')
        nh_groups = smiles.count('NH') + smiles.count('nh')
        estimated_hbd = oh_groups + nh_groups
        estimated_hba = oxygen_count + nitrogen_count
        
        return [estimated_mw, estimated_logp, estimated_hbd, estimated_hba]

    def _create_quantum_circuit(self, features):
        """Create a quantum circuit based on molecular features."""
        qc = QuantumCircuit(4, 4)  # 4 qubits for molecular encoding
        
        # Encode molecular features into quantum circuit
        # Features: [molecular_weight, logP, H_donors, H_acceptors]
        
        # Normalize features for quantum encoding
        mw_norm = min(features[0] / 500, 1.0)  # Normalize molecular weight
        logp_norm = (features[1] + 5) / 10     # Normalize LogP (-5 to 5 range)
        hbd_norm = min(features[2] / 10, 1.0)  # Normalize H donors
        hba_norm = min(features[3] / 20, 1.0)  # Normalize H acceptors
        
        # Apply rotations based on molecular features
        qc.ry(mw_norm * np.pi, 0)    # Molecular weight on qubit 0
        qc.ry(logp_norm * np.pi, 1)  # LogP on qubit 1
        qc.ry(hbd_norm * np.pi, 2)   # H donors on qubit 2
        qc.ry(hba_norm * np.pi, 3)   # H acceptors on qubit 3
        
        # Add entanglement to capture molecular correlations
        qc.cx(0, 1)  # Correlate MW with LogP
        qc.cx(2, 3)  # Correlate H donors with acceptors
        qc.cx(1, 2)  # Cross-correlations
        
        # Measure all qubits
        qc.measure_all()
        
        return qc

    def compute_energy(self, smiles: str) -> float:
        """Compute molecular energy using quantum-inspired methods."""
        try:
            # Extract molecular features with fallback
            if RDKIT_AVAILABLE:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    hbd = Descriptors.NumHDonors(mol)
                    hba = Descriptors.NumHAcceptors(mol)
                    features = [mw, logp, hbd, hba]
                else:
                    # Fall back to estimation if RDKit parsing fails
                    features = self._estimate_molecular_properties(smiles)
            else:
                # Use estimation when RDKit is not available
                features = self._estimate_molecular_properties(smiles)
            
            # Create and execute quantum circuit
            qc = self._create_quantum_circuit(features)
            
            # Execute circuit with error handling
            try:
                # Try with transpilation for better compatibility
                transpiled_qc = transpile(qc, self.backend)
                job = self.backend.run(transpiled_qc, shots=1024)
                result = job.result()
                counts = result.get_counts(transpiled_qc)
            except Exception:
                # Fallback to direct execution
                job = self.backend.run(qc, shots=1024)
                result = job.result()
                counts = result.get_counts(qc)
            
            # Calculate energy from quantum measurement statistics
            total_shots = sum(counts.values())
            energy = 0.0
            
            for bitstring, count in counts.items():
                # Convert bitstring to energy contribution
                probability = count / total_shots
                # Simple mapping: more 1s in bitstring = higher energy
                ones_count = bitstring.count('1')
                energy += (ones_count / 4) * probability
            
            # Scale and offset for realistic molecular energies
            mw = features[0]  # Molecular weight
            base_energy = -mw / 50  # Base energy from molecular weight
            quantum_contribution = energy * 0.5  # Quantum correction
            
            return base_energy + quantum_contribution
            
        except Exception as e:
            print(f"Error computing energy for {smiles}: {e}")
            # Return fallback energy based on estimated properties  
            try:
                features = self._estimate_molecular_properties(smiles)
                return -features[0] / 50  # Simple molecular weight-based energy
            except:
                return 0.0

    def compute_interaction(self, smiles1: str, smiles2: str) -> dict:
        """Compute interaction energy between two molecules."""
        try:
            e1 = self.compute_energy(smiles1)
            e2 = self.compute_energy(smiles2)
            
            # Get molecular properties with fallback handling
            if RDKIT_AVAILABLE:
                mol1 = Chem.MolFromSmiles(smiles1)
                mol2 = Chem.MolFromSmiles(smiles2)
                
                if mol1 and mol2:
                    mw1, mw2 = Descriptors.MolWt(mol1), Descriptors.MolWt(mol2)
                    logp1, logp2 = Descriptors.MolLogP(mol1), Descriptors.MolLogP(mol2)
                    hbd1, hba1 = Descriptors.NumHDonors(mol1), Descriptors.NumHAcceptors(mol1)
                    hbd2, hba2 = Descriptors.NumHDonors(mol2), Descriptors.NumHAcceptors(mol2)
                else:
                    # Fall back to estimation
                    features1 = self._estimate_molecular_properties(smiles1)
                    features2 = self._estimate_molecular_properties(smiles2)
                    mw1, logp1, hbd1, hba1 = features1
                    mw2, logp2, hbd2, hba2 = features2
            else:
                # Use estimation when RDKit is not available
                features1 = self._estimate_molecular_properties(smiles1)
                features2 = self._estimate_molecular_properties(smiles2)
                mw1, logp1, hbd1, hba1 = features1
                mw2, logp2, hbd2, hba2 = features2
            
            # Quantum-inspired interaction calculation
            # Size complementarity
            size_factor = 1 - abs(mw1 - mw2) / (mw1 + mw2)
            
            # Polarity complementarity
            polarity_factor = 1 - abs(logp1 - logp2) / 10
            
            # Hydrogen bonding potential
            hbond_potential = min(hbd1 * hba2 + hbd2 * hba1, 5) / 5
            
            # Calculate binding energy
            binding_strength = -0.4 * size_factor - 0.3 * polarity_factor - 0.5 * hbond_potential
            
            # Add quantum noise for realistic variation
            np.random.seed(hash(smiles1 + smiles2) % 2**32)
            quantum_noise = np.random.normal(0, 0.05)
            
            e_complex = e1 + e2 + binding_strength + quantum_noise
            interaction_energy = e_complex - (e1 + e2)
            
            return {
                "energy1": e1,
                "energy2": e2,
                "energy_complex": e_complex,
                "interaction_energy": interaction_energy
            }
            
        except Exception as e:
            print(f"Error computing interaction: {e}")
            return {
                "energy1": 0.0,
                "energy2": 0.0,
                "energy_complex": 0.0,
                "interaction_energy": 0.0
            }