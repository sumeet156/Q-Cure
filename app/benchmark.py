import time
import numpy as np
from quantum_chemistry import QuantumChemistry
from rdkit import Chem
from rdkit.Chem import Descriptors

def compute_energy_classically(smiles):
    """Classical energy computation using molecular descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 0.0
    
    # Simple classical approximation using molecular weight and properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    
    # Classical energy approximation (simplified)
    classical_energy = -mw / 60 + logp * 0.05
    return classical_energy

def benchmark(smiles):
    """Benchmark quantum vs classical energy computation."""
    print(f"Benchmarking molecular energy computation for: {smiles}")
    print("-" * 50)
    
    # Initialize quantum chemistry
    qc = QuantumChemistry(use_ibmq=False)
    
    # Quantum computation
    start = time.time()
    q_energy = qc.compute_energy(smiles)
    q_time = time.time() - start
    
    # Classical computation
    start = time.time()
    c_energy = compute_energy_classically(smiles)
    c_time = time.time() - start
    
    print(f"Quantum Energy: {q_energy:.4f} Ha (Time: {q_time:.3f}s)")
    print(f"Classical Energy: {c_energy:.4f} Ha (Time: {c_time:.3f}s)")
    
    if q_time > 0:
        speedup = c_time / q_time
        print(f"Speedup: {speedup:.2f}x {'(Quantum faster)' if speedup > 1 else '(Classical faster)'}")
    
    print(f"Energy difference: {abs(q_energy - c_energy):.4f} Ha")
    print()

if __name__ == "__main__":
    # Test on different molecules
    test_molecules = [
        "CCO",              # Ethanol
        "CC(=O)O",          # Acetic acid
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen (simplified)
        "CCC"               # Propane
    ]
    
    for smiles in test_molecules:
        try:
            benchmark(smiles)
        except Exception as e:
            print(f"Error benchmarking {smiles}: {e}")
            print()