#!/usr/bin/env python
"""
qiskit_pipeline.py -- Qiskit stage of the QSP / QSVT pipeline.

Reads exports/qsp_data.json produced by the Wolfram build and:
  * rebuilds the single-qubit QSP gate circuits from the phase angles,
  * cross-checks every circuit unitary against the Wolfram reference,
  * generates OpenQASM 2 and 3,
  * draws circuit diagrams (PNG),
  * simulates with Statevector and the Aer sampler,
  * reconstructs e^{-iHt} from the Qiskit-measured responses,
  * estimates gate resources (native + transpiled + scaling table).

Convention (paper's Wx):
  W(a) = e^{i arccos(a) X} = RX(-2 arccos a)
  S(phi) = e^{i phi Z}     = RZ(-2 phi)
  U = S(phi0) . W . S(phi1) . W ... W . S(phi_d),   P(a) = <0|U|0>
Gates are appended in reverse so the Qiskit unitary equals U.
"""

import json
import os
import numpy as np

from qiskit import QuantumCircuit, transpile, qasm2, qasm3
from qiskit.quantum_info import Operator, Statevector
from qiskit_aer import AerSimulator

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORTS = os.path.join(HERE, "exports")
QASM_DIR = os.path.join(EXPORTS, "qasm")
FIG_DIR = os.path.join(EXPORTS, "figures")
for d in (QASM_DIR, FIG_DIR):
    os.makedirs(d, exist_ok=True)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def c(pair):
    """{re,im} pair -> complex."""
    return complex(pair[0], pair[1])


def cmat(m):
    return np.array([[c(x) for x in row] for row in m], dtype=complex)


def qsp_circuit(phases, a):
    """Single-qubit QSP circuit with U = S(phi0) W S(phi1) ... W S(phi_d)."""
    a = max(-1.0, min(1.0, float(a)))
    theta = np.arccos(a)
    d = len(phases) - 1
    qc = QuantumCircuit(1, name=f"QSP_d{d}")
    # append in reverse: S_d, then (W, S_k) for k = d-1 .. 0
    qc.rz(-2.0 * phases[d], 0)
    for k in range(d - 1, -1, -1):
        qc.rx(-2.0 * theta, 0)
        qc.rz(-2.0 * phases[k], 0)
    return qc


def load_data():
    with open(os.path.join(EXPORTS, "qsp_data.json"), "r") as f:
        return json.load(f)


# ----------------------------------------------------------------------------
# 1. cross-check every checkpoint circuit against the Wolfram reference
# ----------------------------------------------------------------------------
def cross_check(data):
    print("[1] Cross-checking Qiskit circuits vs Wolfram reference unitaries")
    max_u_err = 0.0
    max_p_err = 0.0
    n = 0
    for demo_key in ("demo1_chebyshev", "demo2_cos", "demo2_sin"):
        demo = data[demo_key]
        phases = demo["phases"]
        for cp in demo["checkpoints"]:
            qc = qsp_circuit(phases, cp["a"])
            U_qiskit = Operator(qc).data
            U_ref = cmat(cp["unitary"])
            u_err = np.max(np.abs(U_qiskit - U_ref))
            p_qiskit = U_qiskit[0, 0]
            p_err = abs(p_qiskit - c(cp["P"]))
            max_u_err = max(max_u_err, u_err)
            max_p_err = max(max_p_err, p_err)
            n += 1
    print(f"    checkpoints checked : {n}")
    print(f"    max |U_qiskit - U_wolfram| = {max_u_err:.3e}")
    print(f"    max |P_qiskit - P_wolfram| = {max_p_err:.3e}")
    ok = max_u_err < 1e-9 and max_p_err < 1e-9
    print(f"    cross-check {'PASS' if ok else 'FAIL'}")
    return {"checkpoints": n, "max_unitary_error": max_u_err,
            "max_P_error": max_p_err, "pass": bool(ok)}


# ----------------------------------------------------------------------------
# 2. OpenQASM + diagrams + simulation for the showcase circuits
# ----------------------------------------------------------------------------
def emit_qasm_and_draw(qc, name):
    q2 = os.path.join(QASM_DIR, name + ".qasm")
    with open(q2, "w") as f:
        f.write(qasm2.dumps(qc))
    q3 = os.path.join(QASM_DIR, name + ".qasm3")
    with open(q3, "w") as f:
        f.write(qasm3.dumps(qc))
    png = os.path.join(FIG_DIR, name + "_circuit.png")
    try:
        fig = qc.draw("mpl", fold=40)
        fig.savefig(png, dpi=130, bbox_inches="tight")
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception as e:  # pragma: no cover
        png = None
        print(f"    (diagram skipped for {name}: {e})")
    return q2, q3, png


def showcase_demo1(data):
    print("[2] Demo 1 (QSP -> Chebyshev): showcase circuit, QASM, simulation")
    demo = data["demo1_chebyshev"]
    a = demo["showcase_a"]
    d = demo["degree"]
    qc = qsp_circuit(demo["phases"], a)
    name = f"demo1_chebyshev_d{d}_a{a}"
    q2, q3, png = emit_qasm_and_draw(qc, name)

    sv = Statevector.from_instruction(qc)
    p_amp = sv.data[0]                       # <0|U|0> = P(a) = T_d(a)
    t_d = float(np.polynomial.chebyshev.Chebyshev.basis(d)(a))
    prob0 = abs(p_amp) ** 2

    # Aer shots: probability of measuring |0> = |P|^2
    mc = qc.copy()
    mc.measure_all()
    sim = AerSimulator()
    shots = 20000
    counts = sim.run(transpile(mc, sim), shots=shots).result().get_counts()
    freq0 = counts.get("0", 0) / shots

    print(f"    a={a}, d={d}:  <0|U|0>={p_amp:.6f}   T_{d}({a})={t_d:.6f}")
    print(f"    |P|^2 (exact)={prob0:.5f}   Aer freq(0) over {shots} shots={freq0:.5f}")
    return {
        "a": a, "degree": d, "P_real": p_amp.real, "P_imag": p_amp.imag,
        "chebyshev_T_d": t_d, "abs_err_vs_Td": abs(p_amp.real - t_d),
        "prob0_exact": prob0, "aer_freq0": freq0, "shots": shots,
        "qasm2": os.path.basename(q2), "qasm3": os.path.basename(q3),
        "diagram": os.path.basename(png) if png else None,
    }


def showcase_demo2(data):
    print("[3] Demo 2 (QSP Hamiltonian simulation): per-eigenvalue circuits")
    cos = data["demo2_cos"]
    sin = data["demo2_sin"]
    ham = data["hamiltonian"]
    eta = cos["eta"]
    t = cos["t"]
    eigvals = ham["eigenvalues"]
    eigvecs = cmat(ham["eigenvectors"])      # rows are eigenvectors (Wolfram Eigensystem)
    dim = ham["dimension"]

    # run each eigenvalue's cos/sin circuit on the statevector simulator
    cos_resp, sin_resp = [], []
    for lam in eigvals:
        sc = Statevector.from_instruction(qsp_circuit(cos["phases"], lam))
        ss = Statevector.from_instruction(qsp_circuit(sin["phases"], lam))
        cos_resp.append(sc.data[0].real)
        sin_resp.append(ss.data[0].real)

    # reconstruct e^{-iHt} from the Qiskit-measured responses
    recon = np.zeros((dim, dim), dtype=complex)
    for j, lam in enumerate(eigvals):
        pc = cos_resp[j] / (1 - eta)
        ps = sin_resp[j] / (1 - eta)
        v = eigvecs[j]
        proj = np.outer(v, np.conj(v))
        recon += (pc - 1j * ps) * proj
    exact = cmat(data["reconstruction"]["exact"])
    err = float(np.max(np.abs(recon - exact)))
    overlap = float(abs(np.trace(exact.conj().T @ recon)) / dim)
    print(f"    e^(-iHt) reconstructed from Qiskit outputs:")
    print(f"      max|err| vs exact = {err:.3e}   overlap = {overlap:.10f}")

    # showcase QASM/diagram at the largest eigenvalue
    lam = eigvals[-1]
    qc_cos = qsp_circuit(cos["phases"], lam)
    qc_sin = qsp_circuit(sin["phases"], lam)
    c2, c3, cpng = emit_qasm_and_draw(qc_cos, f"demo2_cos_d{cos['degree']}_eig{lam:.4f}")
    s2, s3, spng = emit_qasm_and_draw(qc_sin, f"demo2_sin_d{sin['degree']}_eig{lam:.4f}")

    return {
        "t": t, "eta": eta,
        "eigenvalues": eigvals,
        "cos_response_qiskit": cos_resp,
        "cos_target": [(1 - eta) * np.cos(t * l) for l in eigvals],
        "sin_response_qiskit": sin_resp,
        "sin_target": [(1 - eta) * np.sin(t * l) for l in eigvals],
        "evolution_max_error": err,
        "evolution_overlap": overlap,
        "cos_qasm2": os.path.basename(c2), "sin_qasm2": os.path.basename(s2),
        "cos_diagram": os.path.basename(cpng) if cpng else None,
        "sin_diagram": os.path.basename(spng) if spng else None,
    }


# ----------------------------------------------------------------------------
# 3. resource estimation
# ----------------------------------------------------------------------------
def resource_row(qc):
    native = dict(qc.count_ops())
    tu = transpile(qc, basis_gates=["u", "cx"], optimization_level=0)
    ibm = transpile(qc, basis_gates=["sx", "rz", "x", "cx"], optimization_level=1)
    opt = transpile(qc, basis_gates=["u", "cx"], optimization_level=3)
    return {
        "num_qubits": qc.num_qubits,
        "native_ops": native,
        "native_depth": qc.depth(),
        "native_size": qc.size(),
        "u_basis_gates": tu.size(), "u_basis_depth": tu.depth(),
        "ibm_basis_ops": dict(ibm.count_ops()), "ibm_basis_depth": ibm.depth(),
        "opt3_u_gates": opt.size(), "opt3_depth": opt.depth(),
    }


def resources(data):
    print("[4] Resource estimation")
    res = {}
    res["demo1_chebyshev_d5"] = resource_row(
        qsp_circuit(data["demo1_chebyshev"]["phases"], data["demo1_chebyshev"]["showcase_a"]))
    res["demo2_cos_d12"] = resource_row(
        qsp_circuit(data["demo2_cos"]["phases"], data["hamiltonian"]["eigenvalues"][-1]))
    res["demo2_sin_d13"] = resource_row(
        qsp_circuit(data["demo2_sin"]["phases"], data["hamiltonian"]["eigenvalues"][-1]))

    # scaling table: degree d -> gate resources (signal calls == degree == queries)
    table = []
    for deg in range(1, 14):
        qc = qsp_circuit([0.1 * (k + 1) for k in range(deg + 1)], 0.6)
        native = qc.count_ops()
        table.append({
            "degree": deg,
            "signal_calls": deg,            # queries to the block encoding W(a)
            "rz": int(native.get("rz", 0)),
            "rx": int(native.get("rx", 0)),
            "total_rotations": qc.size(),
            "depth": qc.depth(),
        })
    res["scaling_table"] = table
    print(f"    demo1 d5 native: {res['demo1_chebyshev_d5']['native_ops']}, depth={res['demo1_chebyshev_d5']['native_depth']}")
    print(f"    demo2 cos d12 native size={res['demo2_cos_d12']['native_size']}, ibm depth={res['demo2_cos_d12']['ibm_basis_depth']}")
    print(f"    scaling table rows: {len(table)} (degree 1..13)")
    return res


# ----------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------
def write_report(cc, d1, d2, res):
    lines = []
    ap = lines.append
    ap("# QSP / QSVT Pipeline -- Qiskit Stage Report\n")
    ap("Circuits are rebuilt from the Wolfram-generated phase angles and verified "
       "against the Wolfram reference unitaries.\n")
    ap("## Cross-tool verification (Wolfram vs Qiskit)\n")
    ap(f"- checkpoints checked: **{cc['checkpoints']}**")
    ap(f"- max |U_qiskit - U_wolfram|: **{cc['max_unitary_error']:.3e}**")
    ap(f"- max |P_qiskit - P_wolfram|: **{cc['max_P_error']:.3e}**")
    ap(f"- result: **{'PASS' if cc['pass'] else 'FAIL'}**\n")

    ap("## Demo 1 - QSP produces Chebyshev polynomials\n")
    ap(f"- showcase a = {d1['a']}, degree d = {d1['degree']}")
    ap(f"- <0|U|0> = {d1['P_real']:.6f} (+{d1['P_imag']:.1e} i)  vs  T_{d1['degree']}({d1['a']}) = {d1['chebyshev_T_d']:.6f}")
    ap(f"- |abs error| = {d1['abs_err_vs_Td']:.3e}")
    ap(f"- Aer sampling: |P|^2 exact = {d1['prob0_exact']:.5f}, freq(0) over {d1['shots']} shots = {d1['aer_freq0']:.5f}")
    ap(f"- OpenQASM: `qasm/{d1['qasm2']}`, `qasm/{d1['qasm3']}`; diagram `figures/{d1['diagram']}`\n")

    ap("## Demo 2 - QSP Hamiltonian simulation e^{-iHt}\n")
    ap(f"- evolution time t = {d2['t']}, subnormalization eta = {d2['eta']}")
    ap(f"- H eigenvalues: {', '.join(f'{v:.4f}' for v in d2['eigenvalues'])}")
    ap("")
    ap("| eigenvalue | cos target | cos (Qiskit) | sin target | sin (Qiskit) |")
    ap("|---|---|---|---|---|")
    for i, lam in enumerate(d2["eigenvalues"]):
        ap(f"| {lam:.4f} | {d2['cos_target'][i]:.5f} | {d2['cos_response_qiskit'][i]:.5f} "
           f"| {d2['sin_target'][i]:.5f} | {d2['sin_response_qiskit'][i]:.5f} |")
    ap("")
    ap(f"- **e^(-iHt) reconstructed from Qiskit outputs: max error = {d2['evolution_max_error']:.3e}, "
       f"overlap = {d2['evolution_overlap']:.10f}**")
    ap(f"- OpenQASM: `qasm/{d2['cos_qasm2']}`, `qasm/{d2['sin_qasm2']}`\n")

    ap("## Resource estimation\n")
    r1 = res["demo1_chebyshev_d5"]
    ap(f"- Demo 1 (d=5): native {r1['native_ops']}, depth {r1['native_depth']}; "
       f"IBM basis depth {r1['ibm_basis_depth']}; optimized single-qubit U-gates {r1['opt3_u_gates']}")
    rc = res["demo2_cos_d12"]
    ap(f"- Demo 2 cos (d=12): native size {rc['native_size']}, depth {rc['native_depth']}; "
       f"IBM basis ops {rc['ibm_basis_ops']}, depth {rc['ibm_basis_depth']}")
    ap("")
    ap("Signal calls (= polynomial degree = queries to the block-encoding of H) dominate cost:")
    ap("")
    ap("| degree | signal calls | Rz | Rx | total rotations | depth |")
    ap("|---|---|---|---|---|---|")
    for row in res["scaling_table"]:
        ap(f"| {row['degree']} | {row['signal_calls']} | {row['rz']} | {row['rx']} "
           f"| {row['total_rotations']} | {row['depth']} |")
    ap("")
    ap("Paper query complexity for Hamiltonian simulation: "
       "Theta( |t| + log(1/eps)/loglog(1/eps) ) queries to the block encoding.")

    with open(os.path.join(EXPORTS, "qiskit_report.md"), "w") as f:
        f.write("\n".join(lines))


def main():
    data = load_data()
    cc = cross_check(data)
    d1 = showcase_demo1(data)
    d2 = showcase_demo2(data)
    res = resources(data)

    with open(os.path.join(EXPORTS, "qiskit_crosscheck.json"), "w") as f:
        json.dump({"cross_check": cc, "demo1": d1,
                   "demo2": {k: v for k, v in d2.items() if k != "eigenvectors"}},
                  f, indent=2)
    with open(os.path.join(EXPORTS, "resources.json"), "w") as f:
        json.dump(res, f, indent=2)
    write_report(cc, d1, d2, res)

    print("\n[done] Qiskit stage complete. Reports in exports/.")
    print(f"  cross-check pass: {cc['pass']}")
    print(f"  demo2 e^(-iHt) overlap: {d2['evolution_overlap']:.10f}")


if __name__ == "__main__":
    main()
