# QSP / QSVT Pipeline -- Qiskit Stage Report

Circuits are rebuilt from the Wolfram-generated phase angles and verified against the Wolfram reference unitaries.

## Cross-tool verification (Wolfram vs Qiskit)

- checkpoints checked: **19**
- max |U_qiskit - U_wolfram|: **6.123e-16**
- max |P_qiskit - P_wolfram|: **3.062e-16**
- result: **PASS**

## Demo 1 - QSP produces Chebyshev polynomials

- showcase a = 0.6, degree d = 5
- <0|U|0> = -0.075840 (+0.0e+00 i)  vs  T_5(0.6) = -0.075840
- |abs error| = 1.665e-16
- Aer sampling: |P|^2 exact = 0.00575, freq(0) over 20000 shots = 0.00635
- OpenQASM: `qasm/demo1_chebyshev_d5_a0.6.qasm`, `qasm/demo1_chebyshev_d5_a0.6.qasm3`; diagram `figures/demo1_chebyshev_d5_a0.6_circuit.png`

## Demo 2 - QSP Hamiltonian simulation e^{-iHt}

- evolution time t = 2.0, subnormalization eta = 0.1
- H eigenvalues: -0.7060, -0.3107, 0.1667, 0.8500

| eigenvalue | cos target | cos (Qiskit) | sin target | sin (Qiskit) |
|---|---|---|---|---|
| -0.7060 | 0.14233 | 0.14233 | -0.88867 | -0.88867 |
| -0.3107 | 0.73180 | 0.73180 | -0.52390 | -0.52390 |
| 0.1667 | 0.85047 | 0.85047 | 0.29445 | 0.29445 |
| 0.8500 | -0.11596 | -0.11596 | 0.89250 | 0.89250 |

- **e^(-iHt) reconstructed from Qiskit outputs: max error = 1.343e-07, overlap = 0.9999999393**
- OpenQASM: `qasm/demo2_cos_d12_eig0.8500.qasm`, `qasm/demo2_sin_d13_eig0.8500.qasm`

## Resource estimation

- Demo 1 (d=5): native {'rz': 6, 'rx': 5}, depth 11; IBM basis depth 5; optimized single-qubit U-gates 1
- Demo 2 cos (d=12): native size 25, depth 25; IBM basis ops {'rz': 3, 'sx': 2}, depth 5

Signal calls (= polynomial degree = queries to the block-encoding of H) dominate cost:

| degree | signal calls | Rz | Rx | total rotations | depth |
|---|---|---|---|---|---|
| 1 | 1 | 2 | 1 | 3 | 3 |
| 2 | 2 | 3 | 2 | 5 | 5 |
| 3 | 3 | 4 | 3 | 7 | 7 |
| 4 | 4 | 5 | 4 | 9 | 9 |
| 5 | 5 | 6 | 5 | 11 | 11 |
| 6 | 6 | 7 | 6 | 13 | 13 |
| 7 | 7 | 8 | 7 | 15 | 15 |
| 8 | 8 | 9 | 8 | 17 | 17 |
| 9 | 9 | 10 | 9 | 19 | 19 |
| 10 | 10 | 11 | 10 | 21 | 21 |
| 11 | 11 | 12 | 11 | 23 | 23 |
| 12 | 12 | 13 | 12 | 25 | 25 |
| 13 | 13 | 14 | 13 | 27 | 27 |

Paper query complexity for Hamiltonian simulation: Theta( |t| + log(1/eps)/loglog(1/eps) ) queries to the block encoding.