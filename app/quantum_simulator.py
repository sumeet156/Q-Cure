from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def simple_quantum_interaction():
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.measure(0, 0)

    simulator = AerSimulator()
    job = simulator.run(qc, shots=1000)
    result = job.result()

    return result.get_counts(qc)
