# References

## Primary papers

1. **J. M. Martyn, Z. M. Rossi, A. K. Tan, I. L. Chuang**,
   *Grand Unification of Quantum Algorithms*,
   PRX Quantum **2**, 040203 (2021). arXiv:2105.02859.
   — The source paper: QSP, QSVT, and how search, simulation, and linear
   systems unify under the quantum eigenvalue/singular-value transform.
   Reproduced here: the QSP theorem (Eq. 4), the zero-phase ⇒ Chebyshev result
   (Eqs. 6–9), block encoding, and the QSVT eigenvalue-transformation picture.

2. **G. H. Low, I. L. Chuang**,
   *Optimal Hamiltonian Simulation by Quantum Signal Processing*,
   Phys. Rev. Lett. **118**, 010501 (2017). arXiv:1606.02685.
   — The QSP Hamiltonian-simulation construction reproduced in Demo 2:
   Jacobi–Anger expansion of `cos`/`sin` (Eqs. 11–12) and the optimal query
   complexity (Eq. 16).

## Methods used

3. **Y. Dong, X. Meng, K. B. Whaley, L. Lin**,
   *Efficient phase-factor evaluation in quantum signal processing*,
   Phys. Rev. A **103**, 042419 (2021). arXiv:2002.11649.
   — The optimization-based phase-angle finder implemented in
   `FindQSPPhases` (symmetric phases, least squares on Chebyshev nodes).

4. **A. Gilyén, Y. Su, G. H. Low, N. Wiebe**,
   *Quantum singular value transformation and beyond*, STOC 2019.
   arXiv:1806.01838. — QSVT formalism underlying the block-encoding /
   eigenvalue-transformation description in `docs/plain_english_derivation.md`.

## Software

5. **Wolfram Quantum Framework** (`Wolfram/QuantumFramework`, v1.6.5), including
   its **Second Quantization** module (`AnnihilationOperator`, `FockState`,
   `$FockSize`, …) used to build the driven-oscillator Hamiltonian.
   Install: `PacletInstall["Wolfram/QuantumFramework"]`;
   the SQ symbols are loaded on demand via
   `PacletSymbol["Wolfram/QuantumFramework", "AnnihilationOperator"]`.
   Docs: <https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/>

6. **Qiskit** 2.5 and **Qiskit Aer** 0.17 — circuit construction, OpenQASM 2/3
   export, statevector and sampling simulation, transpilation and resource
   counting.

## Next step (documented, not implemented here)

7. **Classiq** synthesis — the planned replacement for the Qiskit stage: feed
   the same phase angles / target functions to Classiq's synthesis engine to
   obtain hardware-optimized circuits and a resource comparison. The hand-off
   point is `exports/qsp_data.json` (phases, targets, Hamiltonian).
