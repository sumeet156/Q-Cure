from qiskit.circuit.library import ZZFeatureMap
from qiskit_aer import AerSimulator

def encode_features(features):
    feature_map = ZZFeatureMap(feature_dimension=len(features), reps=1, entanglement='linear')
    qc = feature_map.assign_parameters(features)
    
    # Decompose the circuit to basic gates
    qc = qc.decompose()
    
    # Add save_statevector instruction
    qc.save_statevector()

    simulator = AerSimulator(method='statevector')
    job = simulator.run(qc)
    result = job.result()
    statevector = result.get_statevector()

    return qc, statevector
