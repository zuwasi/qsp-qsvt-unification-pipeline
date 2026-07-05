#!/usr/bin/env wolframscript
(* ============================================================================
   proof_audit.wls -- paper-specific audit of the QSP / QSVT reproduction.
   Symbolic identity checks + validation of the exported artifacts.
   Writes exports/audit_report.txt ; nonzero exit iff a hard check fails.
   ============================================================================ *)

projectDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{projectDir, "src", "QSPUnification.wl"}]];
exportDir = FileNameJoin[{projectDir, "exports"}];

passed = 0; failed = 0; warnings = 0; lines = {};
log[s_] := (AppendTo[lines, s]; Print[s]);
audit[name_, test_, detail_: ""] := If[TrueQ[test],
  (passed++; log["PASS: " <> name]),
  (failed++; log["FAIL: " <> name <> If[detail === "", "", "  -- " <> ToString[detail, InputForm]]])];
warn[name_, detail_] := (warnings++; log["WARN: " <> name <> "  -- " <> ToString[detail, InputForm]]);

log["=== QSP / QSVT Proof Audit ==="];
log["Paper: Grand Unification of Quantum Algorithms (Martyn, Rossi, Tan, Chuang 2021)"];
log[""];

SeedRandom[99];
aRand = RandomReal[{-0.98, 0.98}, 30];

(* ---- A. Convention checks --------------------------------------------------- *)
log["-- A. Operator conventions --"];
audit["W(a)=e^{i arccos(a) X} matches MatrixExp",
  Max@Table[Max@Abs[QSPUnification`QSPSignalOperator[a] -
     MatrixExp[I ArcCos[a] {{0, 1}, {1, 0}}]], {a, aRand}] < 1.*^-12];
audit["S(phi)=e^{i phi Z} matches MatrixExp",
  Max@Table[Max@Abs[QSPUnification`QSPPhaseOperator[p] -
     MatrixExp[I p {{1, 0}, {0, -1}}]], {p, RandomReal[{-Pi, Pi}, 10]}] < 1.*^-12];

(* ---- B. Symbolic Chebyshev derivation (paper eqns) ------------------------- *)
log["-- B. QSP -> Chebyshev symbolic identities --"];
Do[
  With[{cd = QSPUnification`ChebyshevDerivation[d]},
    audit["d=" <> ToString[d] <> ": zero-phase P = T_d", cd["P_equals_Td"]];
    audit["d=" <> ToString[d] <> ": zero-phase Q = U_{d-1}", cd["Q_equals_Udm1"]];
    audit["d=" <> ToString[d] <> ": Pell identity T_d^2+(1-a^2)U_{d-1}^2=1", cd["ChebyshevIdentity"]]],
  {d, {1, 2, 3, 5, 8}}];

(* ---- C. Numeric QSP theorem properties ------------------------------------- *)
log["-- C. QSP theorem: value, unitarity, parity --"];
audit["zero-phase Re P_d(a) = T_d(a) numerically (d<=8)",
  Max@Table[Max@Table[Abs[QSPUnification`QSPRealResponse[ConstantArray[0., d + 1], a] -
     ChebyshevT[d, a]], {a, aRand}], {d, 1, 8}] < 1.*^-11];
audit["QSP unitarity |P|^2+(1-a^2)|Q|^2=1 (random symmetric phases)",
  Module[{ph = RandomReal[{-Pi, Pi}, 4]},
    ph = Join[ph, Reverse[ph]];  (* symmetric, degree 7 *)
    Max@Table[With[{u = QSPUnification`QSPUnitary[ph, a]},
       Abs[Abs[u[[1, 1]]]^2 + Abs[u[[1, 2]]]^2 - 1]], {a, aRand}] < 1.*^-11]];
audit["parity: T_d(-a) = (-1)^d T_d(a) via zero-phase QSP",
  Max@Table[Abs[QSPUnification`QSPRealResponse[ConstantArray[0., d + 1], -0.37] -
     (-1)^d QSPUnification`QSPRealResponse[ConstantArray[0., d + 1], 0.37]], {d, 1, 8}] < 1.*^-11];

(* ---- D. Jacobi-Anger expansion (Hamiltonian simulation kernel) ------------- *)
log["-- D. Jacobi-Anger expansions --"];
With[{t = 2.0, K = 6},
  audit["Jacobi-Anger cos coefficient c_{2k} = 2(-1)^k J_{2k}(t)",
    Module[{ja = QSPUnification`JacobiAngerCos[t, K]},
      Max@Abs[ja["Coefficients"] - Prepend[Table[2 (-1)^k BesselJ[2 k, t], {k, 1, K}], BesselJ[0, t]]] < 1.*^-12]];
  audit["cos(t a) series error < 1e-6 on [-1,1]",
    Max@Table[Abs[QSPUnification`JacobiAngerCos[t, K]["Function"][a] - Cos[t a]], {a, aRand}] < 1.*^-6];
  audit["sin(t a) series error < 1e-6 on [-1,1]",
    Max@Table[Abs[QSPUnification`JacobiAngerSin[t, K]["Function"][a] - Sin[t a]], {a, aRand}] < 1.*^-6]];

(* ---- E. Audit the EXPORTED artifacts (qsp_data.json) ----------------------- *)
log["-- E. Exported artifact validation --"];
jsonPath = FileNameJoin[{exportDir, "qsp_data.json"}];
If[! FileExistsQ[jsonPath],
  warn["qsp_data.json missing", "run build_project.wls first"],
  Module[{data, toC, toM, phc, phs, d1, chk, reconMax, ovl},
    data = Import[jsonPath, "RawJSON"];
    toC[{re_, im_}] := re + I im;
    toM[m_] := Map[toC, m, {2}];

    (* E1: exported zero phases reproduce T_5 *)
    d1 = data["demo1_chebyshev"];
    audit["exported demo1 phases reproduce T_5 (max err < 1e-10)",
      Max[Abs[#["P"][[1]] - #["T_d"]] & /@ d1["checkpoints"]] < 1.*^-10];

    (* E2: exported checkpoints are unitary and match a fresh QSPUnitary build *)
    audit["exported demo1 checkpoint unitaries match recomputation",
      Max[Max@Abs[toM[#["unitary"]] -
         QSPUnification`QSPUnitary[N@d1["phases"], #["a"]]] & /@ d1["checkpoints"]] < 1.*^-10];

    (* E3: exported cos/sin phases reproduce their targets at the eigenvalues *)
    phc = data["demo2_cos"]; phs = data["demo2_sin"];
    audit["exported cos phases hit target at eigenvalues (< 1e-5)",
      Max[Abs[#["response"] - #["target"]] & /@ phc["checkpoints"]] < 1.*^-5];
    audit["exported sin phases hit target at eigenvalues (< 1e-5)",
      Max[Abs[#["response"] - #["target"]] & /@ phs["checkpoints"]] < 1.*^-5];

    (* E4: phase symmetry phi_k = phi_{d-k} *)
    audit["cos phases symmetric",
      Max@Abs[phc["phases"] - Reverse[phc["phases"]]] < 1.*^-6];
    audit["sin phases symmetric",
      Max@Abs[phs["phases"] - Reverse[phs["phases"]]] < 1.*^-6];

    (* E5: reconstruction quality *)
    reconMax = data["reconstruction"]["max_abs_error"];
    ovl = data["reconstruction"]["overlap"];
    audit["e^{-iHt} reconstruction max error < 1e-4", reconMax < 1.*^-4, reconMax];
    audit["e^{-iHt} overlap > 0.9999", ovl > 0.9999, ovl];

    (* E6: exported Hamiltonian is Hermitian with ||H|| <= 1 *)
    Module[{Hm = toM[data["hamiltonian"]["matrix"]]},
      audit["exported H Hermitian", Max@Abs[Hm - ConjugateTranspose[Hm]] < 1.*^-10];
      audit["exported H spectral norm <= 1", Max@Abs@Eigenvalues[Hm] <= 1.0 + 1.*^-9]];
  ]];

(* ---- F. cross-tool (Qiskit) result, if present ---------------------------- *)
log["-- F. Cross-tool (Qiskit) --"];
ccPath = FileNameJoin[{exportDir, "qiskit_crosscheck.json"}];
If[! FileExistsQ[ccPath],
  warn["qiskit_crosscheck.json missing", "run qiskit_pipeline.py to cross-validate"],
  Module[{cc = Import[ccPath, "RawJSON"]},
    audit["Qiskit vs Wolfram unitary agreement < 1e-9",
      cc["cross_check"]["max_unitary_error"] < 1.*^-9, cc["cross_check"]["max_unitary_error"]];
    audit["Qiskit-reconstructed e^{-iHt} overlap > 0.9999",
      cc["demo2"]["evolution_overlap"] > 0.9999, cc["demo2"]["evolution_overlap"]]]];

log[""];
log["=== SUMMARY: " <> ToString[passed] <> " passed, " <> ToString[failed] <>
   " failed, " <> ToString[warnings] <> " warnings ==="];
Export[FileNameJoin[{exportDir, "audit_report.txt"}], StringRiffle[lines, "\n"], "Text"];

If[failed > 0, Exit[1], Exit[0]];
