# QSP / QSVT Pipeline — Grand Unification of Quantum Algorithms

An executable reproduction of the **Quantum Signal Processing (QSP)** core of
J. M. Martyn, Z. M. Rossi, A. K. Tan, I. L. Chuang, *Grand Unification of
Quantum Algorithms*, PRX Quantum **2**, 040203 (2021), together with the QSP
**Hamiltonian-simulation** construction of G. H. Low & I. L. Chuang,
*Optimal Hamiltonian Simulation by Quantum Signal Processing*, PRL **118**,
010501 (2017).

The pipeline runs paper → Mathematica → Qiskit end to end, with a documented
hand-off to Classiq synthesis as the next step.

```
 Paper
   │
   ▼  Wolfram Language (QuantumFramework + SecondQuantization)
 ┌──────────────────────────────────────────────────────────┐
 │ • symbolic: W(a)=e^{iθX},  Wᵈ ⇒ T_d(a),  Pell identity    │
 │ • numerical verification of every identity                 │
 │ • Jacobi–Anger cos/sin (Bessel coefficients)               │
 │ • phase-angle generation (optimization)                    │
 │ • H from second-quantization ops ⇒ diagonalize             │
 │ • reconstruct e^{-iHt}, compare to MatrixExp               │
 └──────────────────────────────────────────────────────────┘
   │  exports/qsp_data.json  (phases, targets, H, unitaries)
   ▼  Qiskit + Aer
 ┌──────────────────────────────────────────────────────────┐
 │ • rebuild single-qubit QSP circuits (Rx/Rz)                │
 │ • cross-check unitary vs Wolfram  (6e-16)                  │
 │ • OpenQASM 2 & 3                                            │
 │ • circuit visualization (PNG)                              │
 │ • statevector + shot-based simulation                      │
 │ • resource estimation (counts / depth / scaling)          │
 └──────────────────────────────────────────────────────────┘
   │
   ▼  (next) Classiq synthesis  — hand-off documented in docs/references.md
```

## Deliverables

| Deliverable | Where |
|---|---|
| **Symbolic derivation** | `src/QSPUnification.wl` (`ChebyshevDerivation`), `docs/plain_english_derivation.md` |
| **Numerical verification** | `exports/validation_report.md`, `exports/validation_results.json`, `exports/audit_report.txt` |
| **Gate circuit** | `qiskit_pipeline.py` (single-qubit QSP), `exports/figures/*_circuit.png` |
| **OpenQASM** | `exports/qasm/*.qasm` (OpenQASM 2) and `*.qasm3` (OpenQASM 3) |
| **Circuit visualization** | `exports/figures/*.png` + Wolfram plots `exports/*.png` |
| **Resource estimation** | `exports/resources.json`, table in `exports/qiskit_report.md` |

## Two demonstrations

* **Demo 1 — QSP ⇒ Chebyshev.** Zero phases make the QSP response exactly the
  Chebyshev polynomial `T_d(a)`. Derived symbolically, verified to machine
  precision, built as a circuit, exported to QASM, and run in Qiskit.
* **Demo 2 — Hamiltonian simulation `e^{-iHt}`.** A driven truncated-Fock
  oscillator (built from the Quantum Framework's second-quantization operators)
  is diagonalized; `cos`/`sin` are approximated by Jacobi–Anger; QSP phase
  angles are generated; the circuits are run in Qiskit and reassembled into
  `e^{-iHt}` (overlap `0.99999994` vs exact).

## File map

```
qsp_pipeline/
├── README.md
├── build_project.wls        # master Wolfram build -> exports/
├── proof_audit.wls          # 35 paper-specific checks -> audit_report.txt
├── qiskit_pipeline.py       # Qiskit stage: circuits, QASM, sim, resources
├── src/
│   └── QSPUnification.wl     # reusable package (QSP core, phases, H, validation)
├── docs/
│   ├── plain_english_derivation.md
│   └── references.md
└── exports/                 # all generated artifacts (see below)
```

## Requirements

* **Wolfram Language 14.x** with paclet `Wolfram/QuantumFramework` (≥1.6.5):
  `PacletInstall["Wolfram/QuantumFramework"]`. Internet is needed once so
  `PacletSymbol[...]` can load the Second Quantization submodule.
* **Python 3** with `qiskit` (≥2.0), `qiskit-aer`, `matplotlib`, `numpy`:
  `pip install qiskit qiskit-aer matplotlib numpy`.

## How to run

```powershell
# 1. Wolfram build: symbolic derivation, verification, phases, H, exports
wolframscript -file build_project.wls

# 2. Qiskit stage: circuits, cross-check, OpenQASM, simulation, resources
python qiskit_pipeline.py

# 3. Paper-specific audit — runs last so it validates both the Wolfram
#    artifacts and the Qiskit cross-check (exits nonzero on any hard failure)
wolframscript -file proof_audit.wls
```

Run order matters: `build_project.wls` writes `exports/qsp_data.json` (consumed
by the other two), and `qiskit_pipeline.py` writes `exports/qiskit_crosscheck.json`
(audited in step 3). The audit still runs standalone — it just emits warnings for
any artifact not yet generated.

## Outputs (`exports/`)

```
qsp_data.json              phases, targets, H, expected unitaries (Wolfram → Qiskit)
validation_results.json    machine-readable validation
validation_report.md       human-readable validation table
audit_report.txt           35 proof-audit checks (all PASS)
qiskit_report.md           cross-check, simulation, resource tables
qiskit_crosscheck.json     Wolfram-vs-Qiskit agreement + reconstruction
resources.json             gate counts / depth / scaling table
chebyshev_verification.csv, phase_angles.csv, eigenvalue_transform.csv
chebyshev_qsp.png, cos_qsp.png, sin_qsp.png,
eigenvalue_transform.png, reconstruction.png       (Wolfram figures)
figures/*_circuit.png                              (Qiskit circuit diagrams)
qasm/*.qasm, qasm/*.qasm3                           (OpenQASM 2 & 3)
```

## Results at a glance

| Check | Result |
|---|---|
| Chebyshev derivation `T_d`, `d=1..8` (symbolic) | exact |
| Pell identity `T_d²+(1-a²)U_{d-1}²=1` | exact |
| zero-phase QSP vs `T_d` (numeric) | `< 10⁻¹¹` |
| Jacobi–Anger cos/sin truncation (`t=2,K=6`) | `< 10⁻⁶` |
| phase generation (cos deg 12 / sin deg 13) | `≈ 3·10⁻⁷ / 7·10⁻⁸` |
| `e^{-iHt}` reconstruction / overlap | `1.3·10⁻⁷` / `0.99999994` |
| **Wolfram ↔ Qiskit unitary agreement** | **`6·10⁻¹⁶`** |
| proof audit | **35 passed, 0 failed** |

## Conventions

Paper's Wx convention throughout:
`W(a)=e^{i·arccos(a)·X}`, `S(φ)=e^{i·φ·Z}`,
`U_φ = S(φ₀)·∏ₖ W(a)·S(φₖ)`, `P(a)=⟨0|U_φ|0⟩`.
Qiskit mapping: `W(a)=RX(-2·arccos a)`, `S(φ)=RZ(-2·φ)`.
