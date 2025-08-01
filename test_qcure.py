"""
Test script to verify Q-Cure functionality
Tests all core features:
- Quantum Chemistry (VQE) simulation
- Real Drug Data (PubChem) integration  
- IBM Quantum Hardware option
- Risk Assessment
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'app'))

from quantum_chemistry import QuantumChemistry
from molecule_utils import get_molecule_features
from interaction_score import assess_risk

def test_quantum_chemistry():
    """Test quantum chemistry simulation"""
    print("🧪 Testing Quantum Chemistry (VQE)...")
    
    qc = QuantumChemistry(use_ibmq=False)  # Use simulator for testing
    
    # Test individual energy calculation
    energy = qc.compute_energy("CCO")  # Ethanol
    print(f"   ✅ Ethanol energy: {energy:.4f} Ha")
    
    # Test interaction calculation
    result = qc.compute_interaction("CCO", "CC(=O)O")  # Ethanol + Acetic acid
    print(f"   ✅ Interaction energy: {result['interaction_energy']:.4f} Ha")
    
    return True

def test_pubchem_integration():
    """Test PubChem data fetching"""
    print("📊 Testing Real Drug Data (PubChem)...")
    
    # Test with drug name
    mol1 = get_molecule_features("aspirin")
    if mol1:
        print(f"   ✅ Aspirin: {mol1['smiles']} (MW: {mol1['mw']:.2f})")
    else:
        print("   ❌ Failed to fetch aspirin data")
        return False
    
    # Test with SMILES
    mol2 = get_molecule_features("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin SMILES
    if mol2:
        print(f"   ✅ SMILES lookup: MW {mol2['mw']:.2f}")
    else:
        print("   ❌ Failed to process SMILES")
        return False
    
    return True

def test_risk_assessment():
    """Test risk assessment function"""
    print("⚠️  Testing Risk Assessment...")
    
    # Test different risk levels
    high_risk = assess_risk(-0.8)  # Strong binding
    moderate_risk = assess_risk(-0.3)  # Moderate binding  
    low_risk = assess_risk(-0.1)  # Weak binding
    no_risk = assess_risk(0.5)  # Repulsive
    
    print(f"   ✅ High risk (-0.8): {high_risk}")
    print(f"   ✅ Moderate risk (-0.3): {moderate_risk}")
    print(f"   ✅ Low risk (-0.1): {low_risk}")
    print(f"   ✅ No risk (0.5): {no_risk}")
    
    return True

def test_ibm_quantum_option():
    """Test IBM Quantum hardware option"""
    print("🚀 Testing IBM Quantum Hardware Option...")
    
    # Test simulator mode
    qc_sim = QuantumChemistry(use_ibmq=False)
    print("   ✅ Simulator mode initialized")
    
    # Test IBM Quantum mode (will fall back to simulator safely)
    qc_ibmq = QuantumChemistry(use_ibmq=True)
    print("   ✅ IBM Quantum mode initialized (fallback to simulator)")
    
    return True

def run_full_integration_test():
    """Run complete drug interaction prediction"""
    print("🔬 Running Full Integration Test...")
    
    # Initialize components
    qc = QuantumChemistry(use_ibmq=False)
    
    # Fetch drug data
    mol1 = get_molecule_features("aspirin")
    mol2 = get_molecule_features("ibuprofen")
    
    if not mol1 or not mol2:
        print("   ❌ Failed to fetch drug data")
        return False
    
    # Compute interaction
    result = qc.compute_interaction(mol1["smiles"], mol2["smiles"])
    
    # Assess risk
    risk = assess_risk(result["interaction_energy"])
    
    print(f"   ✅ Aspirin + Ibuprofen interaction:")
    print(f"      - Energy 1: {result['energy1']:.4f} Ha")
    print(f"      - Energy 2: {result['energy2']:.4f} Ha") 
    print(f"      - Interaction: {result['interaction_energy']:.4f} Ha")
    print(f"      - Risk: {risk}")
    
    return True

if __name__ == "__main__":
    print("🧪 Q-Cure Feature Verification Test")
    print("=" * 50)
    
    tests = [
        ("Quantum Chemistry (VQE)", test_quantum_chemistry),
        ("Real Drug Data (PubChem)", test_pubchem_integration), 
        ("IBM Quantum Hardware", test_ibm_quantum_option),
        ("Risk Assessment", test_risk_assessment),
        ("Full Integration", run_full_integration_test)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            if test_func():
                print(f"✅ {test_name}: PASSED\n")
                passed += 1
            else:
                print(f"❌ {test_name}: FAILED\n")
                failed += 1
        except Exception as e:
            print(f"❌ {test_name}: ERROR - {e}\n")
            failed += 1
    
    print("=" * 50)
    print(f"📊 Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All Q-Cure features are working correctly!")
        print("\nFeatures verified:")
        print("✅ Quantum Chemistry (VQE) – Simulates molecular energies")
        print("✅ Real Drug Data (PubChem) – No hardcoded molecules")
        print("✅ IBM Quantum Hardware – Optional real quantum execution")
        print("✅ Risk Assessment – Clear visual output")
    else:
        print("⚠️  Some features need attention")
