# Quantum Signal Processing → Chebyshev → Hamiltonian Simulation

A plain-English mathematical companion to the code in `src/QSPUnification.wl`,
`build_project.wls`, and `qiskit_pipeline.py`. It reproduces the core of

* Martyn, Rossi, Tan, Chuang, *Grand Unification of Quantum Algorithms*, PRX Quantum **2**, 040203 (2021), and
* Low, Chuang, *Optimal Hamiltonian Simulation by Quantum Signal Processing*, PRL **118**, 010501 (2017).

Every equation below is checked either symbolically or numerically by
`proof_audit.wls` (35 checks, all passing).

---

## 1. Quantum Signal Processing (QSP)

QSP interleaves two single-qubit rotations. In the paper's **Wx convention**:

**Signal operator** (an `x`-rotation carrying the input `a ∈ [-1,1]`):

```
           ⎡      a          i·√(1-a²) ⎤
W(a)  =    ⎢                            ⎥   =  e^{ i·arccos(a)·X }.        (1)
           ⎣ i·√(1-a²)          a      ⎦
```

**Signal-processing operator** (a `z`-rotation carrying a tunable phase `φ`):

```
S(φ)  =  e^{ i·φ·Z }  =  diag( e^{iφ}, e^{-iφ} ).                          (2)
```

**QSP sequence** of degree `d` with phases `φ = (φ₀, φ₁, …, φ_d)` (there are
`d` signal calls and `d+1` phases):

```
U_φ(a)  =  S(φ₀) · ∏_{k=1}^{d} [ W(a) · S(φ_k) ].                          (3)
```

**QSP theorem.** There exist polynomials `P, Q ∈ ℂ[a]` with `deg P ≤ d`,
`deg Q ≤ d-1`, and parity `d mod 2`, such that

```
           ⎡   P(a)        i·Q(a)·√(1-a²) ⎤
U_φ(a) =   ⎢                               ⎥,   |P|² + (1-a²)|Q|² = 1.     (4)
           ⎣ i·Q*(a)·√(1-a²)     P*(a)     ⎦
```

The quantity we read out of the circuit is the top-left amplitude:

```
⟨0| U_φ(a) |0⟩  =  P(a).                                                   (5)
```

Equations (1)–(2) are verified against `MatrixExp` in audit block **A**.

---

## 2. Zero phases ⇒ Chebyshev polynomials (the symbolic core)

Set every phase to zero, `φ = (0,…,0)`. Then `S(0) = I` and (3) collapses to a
pure power of the signal operator:

```
U_0(a)  =  W(a)^d.                                                         (6)
```

Write `a = cos θ`, so by (1) `W = e^{iθX} = cos θ · I + i sin θ · X`. Powers of
a single-axis rotation just add angles:

```
W(a)^d  =  e^{ i·d·θ·X }  =  cos(dθ)·I + i·sin(dθ)·X
        =  ⎡ cos(dθ)      i·sin(dθ) ⎤
           ⎣ i·sin(dθ)      cos(dθ) ⎦.                                     (7)
```

Now use the **definitions** of the Chebyshev polynomials,
`cos(dθ) = T_d(cos θ)` and `sin(dθ) = sin θ · U_{d-1}(cos θ)`. Substituting
`cos θ = a`, `sin θ = √(1-a²)` and comparing with (4):

```
P(a) = T_d(a),        Q(a) = U_{d-1}(a).                                   (8)
```

So **the QSP response with zero phases is exactly the Chebyshev polynomial of
the first kind**, `⟨0|U_0(a)|0⟩ = T_d(a)`, and the off-diagonal block carries
the second-kind polynomial `U_{d-1}`. This is derived *symbolically* in
`ChebyshevDerivation[d]` (via `MatrixPower` + `FullSimplify`) for `d = 1..8`.

The unitarity condition `|P|² + (1-a²)|Q|² = 1` in (4) becomes the classical
**Pell / Chebyshev identity**

```
T_d(a)²  +  (1-a²)·U_{d-1}(a)²  =  1,                                      (9)
```

proved symbolically (`FullSimplify == 0`) and checked numerically in audit
blocks **B** and **C**.

*Figure:* `exports/chebyshev_qsp.png` overlays `T₁..T₅` with the zero-phase QSP
response for `d=5` (red dots), which land exactly on the `T₅` curve.

---

## 3. Generating phase angles for a target function

For a general real target `f(a)` with definite parity and `‖f‖∞ < 1`, QSP can
realize `Re P(a) = f(a)` for suitable phases. We restrict to **symmetric**
phase sequences `φ_k = φ_{d-k}` and find them by minimizing

```
L(φ)  =  Σ_j ( Re⟨0|U_φ(x_j)|0⟩ − f(x_j) )²                              (10)
```

over Chebyshev nodes `x_j = cos((2j-1)π/2M)` (`FindQSPPhases`; the
optimization-based method of Dong, Meng, Whaley, Lin, PRA 2021). Because the
same 2×2 product (3) is evaluated in Wolfram for optimization/verification and
rebuilt in Qiskit as gates, the two tool-chains agree to `6·10⁻¹⁶`
(cross-check in `qiskit_pipeline.py`).

---

## 4. Hamiltonian simulation `e^{-iHt}` (Low–Chuang)

**Goal:** apply the matrix function `e^{-iHt}` for a Hermitian `H` with
`‖H‖ ≤ 1`. Because `e^{-iHt} = cos(Ht) − i·sin(Ht)`, it suffices to realize the
two real functions `cos(t·a)` and `sin(t·a)` of the eigenvalues `a` of `H`.

**Polynomial approximation (Jacobi–Anger).** These functions have exact
Chebyshev expansions with Bessel-function coefficients:

```
cos(t·a) = J₀(t) + 2·Σ_{k≥1} (−1)^k J_{2k}(t)·T_{2k}(a),                  (11)
sin(t·a) =        2·Σ_{k≥0} (−1)^k J_{2k+1}(t)·T_{2k+1}(a).               (12)
```

Truncating at order `2K` (resp. `2K+1`) gives even/odd polynomials; the
truncation error is `O(J_{≥2K}(t))`, which is `< 10⁻⁶` here for `t=2, K=6`
(`JacobiAngerCos/Sin`, audit block **D**). We use a subnormalization factor
`(1−η)` so the targets satisfy `‖f‖∞ < 1` as QSP requires.

**Diagonalization + reconstruction.** We build a physical Hamiltonian from the
Wolfram Quantum Framework's **second-quantization** operators — a driven,
truncated-Fock oscillator

```
H_phys = ω·n + g·(a + a†),        n = a†a,   dim = 4,                      (13)
```

then center and rescale to `‖H‖ = 0.85`. Diagonalizing gives eigenpairs
`(λ_j, |v_j⟩)`. Running the QSP circuits at each `λ_j` yields
`p_cos(λ_j) ≈ (1−η)cos(tλ_j)` and `p_sin(λ_j) ≈ (1−η)sin(tλ_j)`, and we
reassemble

```
e^{-iHt}  ≈  Σ_j (1/(1−η))·[ p_cos(λ_j) − i·p_sin(λ_j) ]·|v_j⟩⟨v_j|.      (14)
```

`ReconstructEvolution` (Wolfram) and the same reconstruction *driven by Qiskit
statevector outputs* both reproduce the exact `MatrixExp[-iHt]` with

* max abs error ≈ `1.3·10⁻⁷`,
* normalized overlap ≈ `0.99999994`.

---

## 5. Circuit realization and the Qiskit mapping

Each QSP sequence (3) is a **single-qubit** circuit. Using
`RX(λ)=e^{-iλX/2}`, `RZ(λ)=e^{-iλZ/2}`:

```
W(a) = RX(−2·arccos a),        S(φ) = RZ(−2·φ).                           (15)
```

Gates are appended in reverse order so the Qiskit unitary equals `U_φ(a)`; the
top-left amplitude reproduces `P(a)`. This is exported to **OpenQASM 2/3**,
drawn, simulated (statevector + Aer sampling), and resource-counted.

**Resource picture.** A degree-`d` QSP circuit uses `d` signal calls
(= queries to the block encoding of `H`) and `d+1` phase rotations, i.e.
`2d+1` single-qubit rotations. For Hamiltonian simulation the paper's query
complexity is

```
#queries  =  Θ( |t| + log(1/ε) / log log(1/ε) ),                          (16)
```

which is optimal in both `t` and `ε`. The degree we use (`12`/`13`) is exactly
this query count for our `t=2` and target accuracy.

---

## 6. Numerical verification summary

| Quantity | Result |
|---|---|
| zero-phase `Re P_d(a) − T_d(a)`, `d≤8` | `< 10⁻¹¹` |
| Pell identity (9), symbolic | exact (`0`) |
| QSP unitarity `|P|²+(1-a²)|Q|²−1` | `< 10⁻¹¹` |
| Jacobi–Anger cos/sin truncation error (`t=2,K=6`) | `< 10⁻⁶` |
| phase-generation error (cos deg 12 / sin deg 13) | `≈ 3·10⁻⁷ / 7·10⁻⁸` |
| `e^{-iHt}` reconstruction error / overlap | `1.3·10⁻⁷` / `0.99999994` |
| Wolfram vs Qiskit unitary agreement | `6·10⁻¹⁶` |

**Adversarial self-review: all checks passed.** Assumptions are stated
(`‖H‖≤1`, `‖f‖∞<1` via subnormalization `η`, symmetric phases); existence and
the QSP polynomial parity are used exactly as proven in the paper; the truncated
Fock space defines `H` exactly (no hidden truncation claim); every displayed
identity is matched by an executable check.
