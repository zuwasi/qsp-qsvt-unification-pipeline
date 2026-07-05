#!/usr/bin/env wolframscript
(* ============================================================================
   build_project.wls  --  master Wolfram build for the QSP / QSVT pipeline
   Reproduces "Grand Unification of Quantum Algorithms" (Martyn et al. 2021)
   and the QSP Hamiltonian-simulation construction (Low & Chuang 2017).

   Produces (in exports/):
     qsp_data.json           phases, targets, H, expected unitaries  (-> Qiskit)
     validation_results.json  machine-readable validation
     validation_report.md     human-readable validation
     chebyshev_verification.csv, phase_angles.csv, eigenvalue_transform.csv
     chebyshev_qsp.png, cos_qsp.png, sin_qsp.png,
     eigenvalue_transform.png, reconstruction.png
   ============================================================================ *)

projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];
Get[FileNameJoin[{projectDir, "src", "QSPUnification.wl"}]];

exportDir = FileNameJoin[{projectDir, "exports"}];
If[! DirectoryQ[exportDir], CreateDirectory[exportDir]];
If[! DirectoryQ[FileNameJoin[{projectDir, "docs"}]],
   CreateDirectory[FileNameJoin[{projectDir, "docs"}]]];

cx = QSPUnification`ComplexToPairs;   (* complex -> {re,im}, recurses over lists *)
ck[x_] := N[Chop[x, 1.*^-14]];
Print["=================================================================="];
Print[" QSP / QSVT pipeline build"];
Print["=================================================================="];

(* ---------------------------------------------------------------------------
   DEMO 1 : QSP  ->  Chebyshev polynomials  (exact, single qubit)
   --------------------------------------------------------------------------- *)
Print["\n[Demo 1] QSP -> Chebyshev polynomials (symbolic + numeric)"];
d1 = 5;
d1Phases = ConstantArray[0., d1 + 1];
symDeriv = Table[QSPUnification`ChebyshevDerivation[d], {d, 1, 6}];
Do[Print["  d=", sd["Degree"], "  P=T_d:", sd["P_equals_Td"],
   "  Q=U_{d-1}:", sd["Q_equals_Udm1"], "  Pell id:", sd["ChebyshevIdentity"]],
 {sd, symDeriv}];

d1Grid = Range[-1., 1., 0.2];
d1Checkpoints = Table[
  With[{u = QSPUnification`QSPUnitary[d1Phases, a]},
    <|"a" -> a, "unitary" -> cx[ck[u]], "P" -> cx[ck[u[[1, 1]]]],
      "T_d" -> N[ChebyshevT[d1, a]]|>],
  {a, d1Grid}];
d1MaxErr = Max@Table[Abs[cp["P"][[1]] - cp["T_d"]], {cp, d1Checkpoints}];
Print["  max |Re P_5(a) - T_5(a)| over grid = ", d1MaxErr];

(* ---------------------------------------------------------------------------
   DEMO 2 : QSP Hamiltonian simulation  e^{-iHt}
   --------------------------------------------------------------------------- *)
Print["\n[Demo 2] QSP Hamiltonian simulation"];
tEvol = 2.0; eta = 0.1; Kexp = 6;
dCos = 2 Kexp; dSin = 2 Kexp + 1;

Print["  building H from second-quantization operators ..."];
H = QSPUnification`BuildHamiltonian[];
eigs = H["Eigenvalues"];
Print["  H dim=", H["Dimension"], "  eigenvalues=", eigs];

Print["  generating cos phases (degree ", dCos, ") ..."];
cosTarget = Function[x, (1 - eta) Cos[tEvol x]];
cosD = QSPUnification`FindQSPPhases[cosTarget, dCos, "Restarts" -> 6];
Print["    cos phase max error = ", cosD["MaxError"]];
Print["  generating sin phases (degree ", dSin, ") ..."];
sinTarget = Function[x, (1 - eta) Sin[tEvol x]];
sinD = QSPUnification`FindQSPPhases[sinTarget, dSin, "Restarts" -> 6];
Print["    sin phase max error = ", sinD["MaxError"]];

mkCheckpoint[phases_, tf_][a_] := With[{u = QSPUnification`QSPUnitary[phases, a]},
  <|"a" -> N[a], "unitary" -> cx[ck[u]], "P" -> cx[ck[u[[1, 1]]]],
    "response" -> N[Re[u[[1, 1]]]], "target" -> N[tf[a]]|>];

cosCheckpoints = mkCheckpoint[cosD["Phases"], cosTarget] /@ eigs;
sinCheckpoints = mkCheckpoint[sinD["Phases"], sinTarget] /@ eigs;

Print["  reconstructing e^{-iHt} ..."];
recon = QSPUnification`ReconstructEvolution[H["Matrix"], tEvol, cosD["Phases"], sinD["Phases"], eta];
Print["    reconstruction max|err| = ", recon["MaxAbsError"],
   "   overlap = ", recon["NormalizedOverlap"]];

(* ---------------------------------------------------------------------------
   Assemble JSON payload for the Qiskit stage
   --------------------------------------------------------------------------- *)
Print["\n[Export] qsp_data.json"];
payload = <|
  "meta" -> <|
    "paper" -> "Grand Unification of Quantum Algorithms (Martyn, Rossi, Tan, Chuang 2021)",
    "companion" -> "Optimal Hamiltonian Simulation by QSP (Low, Chuang 2017)",
    "convention" -> "Wx: W(a)=e^{i arccos(a) X}, S(phi)=e^{i phi Z}, U=S(phi0) Prod W(a) S(phi_k)",
    "qiskit_mapping" -> "W(a)=RX(-2 arccos a),  S(phi)=RZ(-2 phi),  P(a)=<0|U|0>",
    "generated" -> DateString["ISODateTime"]
  |>,
  "demo1_chebyshev" -> <|
    "degree" -> d1, "phases" -> N[d1Phases],
    "checkpoints" -> d1Checkpoints,
    "showcase_a" -> 0.6, "max_error" -> d1MaxErr
  |>,
  "demo2_cos" -> <|
    "t" -> tEvol, "eta" -> eta, "degree" -> dCos, "parity" -> "even",
    "phases" -> N[cosD["Phases"]], "max_error" -> cosD["MaxError"],
    "checkpoints" -> cosCheckpoints
  |>,
  "demo2_sin" -> <|
    "t" -> tEvol, "eta" -> eta, "degree" -> dSin, "parity" -> "odd",
    "phases" -> N[sinD["Phases"]], "max_error" -> sinD["MaxError"],
    "checkpoints" -> sinCheckpoints
  |>,
  "hamiltonian" -> <|
    "dimension" -> H["Dimension"], "t" -> tEvol, "params" -> H["Params"],
    "matrix" -> cx[ck[H["Matrix"]]],
    "eigenvalues" -> N[eigs],
    "eigenvectors" -> cx[ck[H["Eigenvectors"]]],
    "annihilation" -> cx[ck[H["AnnihilationMatrix"]]],
    "number" -> cx[ck[H["NumberMatrix"]]]
  |>,
  "reconstruction" -> <|
    "exact" -> cx[ck[recon["Exact"]]],
    "reconstructed" -> cx[ck[recon["Reconstructed"]]],
    "max_abs_error" -> recon["MaxAbsError"],
    "overlap" -> recon["NormalizedOverlap"]
  |>
|>;
Export[FileNameJoin[{exportDir, "qsp_data.json"}], payload, "JSON"];

(* ---------------------------------------------------------------------------
   Validation suite -> JSON + Markdown
   --------------------------------------------------------------------------- *)
Print["\n[Export] validation results"];
val = QSPUnification`RunValidationSuite[];
chebOK = AllTrue[Values[val["ChebyshevDerivation"]],
  #["P_equals_Td"] && #["Q_equals_Udm1"] && #["ChebyshevIdentity"] &];
valFlat = <|
  "chebyshev_symbolic_all_pass" -> chebOK,
  "zero_phase_chebyshev_max_error" -> Max[Values[val["ZeroPhaseChebyshev"]]],
  "qsp_unitarity_max_error" -> val["QSPUnitarity"],
  "jacobi_anger_cos_error" -> val["JacobiAngerCosError"],
  "jacobi_anger_sin_error" -> val["JacobiAngerSinError"],
  "phase_gen_cos_max_error" -> val["PhaseGenCosMaxError"],
  "phase_gen_sin_max_error" -> val["PhaseGenSinMaxError"],
  "hamiltonian_sim_max_error" -> val["HamiltonianSimError"],
  "hamiltonian_sim_overlap" -> val["HamiltonianSimOverlap"]
|>;
Export[FileNameJoin[{exportDir, "validation_results.json"}], valFlat, "JSON"];

vlines = {
  "# QSP / QSVT Pipeline -- Validation Report", "",
  "Paper: Grand Unification of Quantum Algorithms (Martyn, Rossi, Tan, Chuang, 2021).", "",
  "| Check | Result |", "|---|---|",
  "| Chebyshev symbolic derivation (d=1..6) all pass | " <> ToString[chebOK] <> " |",
  "| Zero-phase QSP vs T_d, max error | " <> ToString[valFlat["zero_phase_chebyshev_max_error"]] <> " |",
  "| QSP unitarity |P|^2+(1-a^2)|Q|^2=1, max error | " <> ToString[valFlat["qsp_unitarity_max_error"]] <> " |",
  "| Jacobi-Anger cos truncation error | " <> ToString[valFlat["jacobi_anger_cos_error"]] <> " |",
  "| Jacobi-Anger sin truncation error | " <> ToString[valFlat["jacobi_anger_sin_error"]] <> " |",
  "| Phase generation, cos target max error | " <> ToString[valFlat["phase_gen_cos_max_error"]] <> " |",
  "| Phase generation, sin target max error | " <> ToString[valFlat["phase_gen_sin_max_error"]] <> " |",
  "| e^{-iHt} reconstruction max error | " <> ToString[valFlat["hamiltonian_sim_max_error"]] <> " |",
  "| e^{-iHt} normalized overlap | " <> ToString[valFlat["hamiltonian_sim_overlap"]] <> " |"
};
Export[FileNameJoin[{exportDir, "validation_report.md"}], StringRiffle[vlines, "\n"], "Text"];

(* ---------------------------------------------------------------------------
   CSV exports
   --------------------------------------------------------------------------- *)
Print["[Export] CSV tables"];
chebCsv = Prepend[
  Table[{sd["Degree"], sd["P_equals_Td"], sd["Q_equals_Udm1"], sd["ChebyshevIdentity"],
     Max@Table[Abs[QSPUnification`QSPRealResponse[ConstantArray[0., sd["Degree"] + 1], a] -
        ChebyshevT[sd["Degree"], a]], {a, Range[-1, 1, 0.1]}]},
    {sd, symDeriv}],
  {"degree", "P=T_d", "Q=U_{d-1}", "PellIdentity", "max_numeric_error"}];
Export[FileNameJoin[{exportDir, "chebyshev_verification.csv"}], chebCsv];

phaseCsv = Join[
  {{"demo", "k", "phase"}},
  Table[{"chebyshev_d5", k - 1, d1Phases[[k]]}, {k, 1, Length[d1Phases]}],
  Table[{"cos_d" <> ToString[dCos], k - 1, cosD["Phases"][[k]]}, {k, 1, Length[cosD["Phases"]]}],
  Table[{"sin_d" <> ToString[dSin], k - 1, sinD["Phases"][[k]]}, {k, 1, Length[sinD["Phases"]]}]];
Export[FileNameJoin[{exportDir, "phase_angles.csv"}], phaseCsv];

eigCsv = Prepend[
  Table[{eigs[[j]],
     N[(1 - eta) Cos[tEvol eigs[[j]]]], cosCheckpoints[[j]]["response"],
     N[(1 - eta) Sin[tEvol eigs[[j]]]], sinCheckpoints[[j]]["response"]},
    {j, 1, Length[eigs]}],
  {"eigenvalue", "cos_target", "cos_qsp", "sin_target", "sin_qsp"}];
Export[FileNameJoin[{exportDir, "eigenvalue_transform.csv"}], eigCsv];

(* ---------------------------------------------------------------------------
   Figures
   --------------------------------------------------------------------------- *)
Print["[Export] figures (PNG)"];
fineGrid = Range[-1., 1., 0.01];

fig1 = Plot[Evaluate@Table[ChebyshevT[d, a], {d, 1, 5}], {a, -1, 1},
  PlotLegends -> Table["T_" <> ToString[d], {d, 1, 5}],
  Epilog -> {PointSize[0.012], Red,
    Point[Table[{a, QSPUnification`QSPRealResponse[ConstantArray[0., 6], a]}, {a, Range[-1, 1, 0.1]}]]},
  PlotLabel -> "Demo 1: T_d(a) (lines) vs zero-phase QSP response for d=5 (red)",
  AxesLabel -> {"a", "P(a)"}, ImageSize -> 640];
Export[FileNameJoin[{exportDir, "chebyshev_qsp.png"}], fig1];

fig2 = Show[
  Plot[(1 - eta) Cos[tEvol a], {a, -1, 1}, PlotStyle -> {Blue, Thick},
    PlotLabel -> "Demo 2: (1-eta) cos(t a) target vs QSP degree-" <> ToString[dCos] <> " response",
    AxesLabel -> {"a", "response"}],
  ListPlot[Table[{a, QSPUnification`QSPRealResponse[cosD["Phases"], a]}, {a, Range[-1, 1, 0.05]}],
    PlotStyle -> {Red, PointSize[0.01]}], ImageSize -> 640];
Export[FileNameJoin[{exportDir, "cos_qsp.png"}], fig2];

fig3 = Show[
  Plot[(1 - eta) Sin[tEvol a], {a, -1, 1}, PlotStyle -> {Darker[Green], Thick},
    PlotLabel -> "Demo 2: (1-eta) sin(t a) target vs QSP degree-" <> ToString[dSin] <> " response",
    AxesLabel -> {"a", "response"}],
  ListPlot[Table[{a, QSPUnification`QSPRealResponse[sinD["Phases"], a]}, {a, Range[-1, 1, 0.05]}],
    PlotStyle -> {Red, PointSize[0.01]}], ImageSize -> 640];
Export[FileNameJoin[{exportDir, "sin_qsp.png"}], fig3];

fig4 = Show[
  Plot[{(1 - eta) Cos[tEvol a], (1 - eta) Sin[tEvol a]}, {a, -1, 1},
    PlotStyle -> {Blue, Darker[Green]},
    PlotLegends -> {"(1-eta)cos(t a)", "(1-eta)sin(t a)"},
    PlotLabel -> "Eigenvalue transform: markers at H eigenvalues",
    AxesLabel -> {"eigenvalue a", "f(a)"}],
  ListPlot[{Table[{eigs[[j]], cosCheckpoints[[j]]["response"]}, {j, Length[eigs]}],
            Table[{eigs[[j]], sinCheckpoints[[j]]["response"]}, {j, Length[eigs]}]},
    PlotStyle -> {{Red, PointSize[0.02]}, {Orange, PointSize[0.02]}}], ImageSize -> 640];
Export[FileNameJoin[{exportDir, "eigenvalue_transform.png"}], fig4];

fig5 = MatrixPlot[Abs[recon["Reconstructed"] - recon["Exact"]],
  PlotLabel -> "| QSP reconstruction  -  exact e^{-iHt} |  (abs error)",
  ColorFunction -> "TemperatureMap", ImageSize -> 480];
Export[FileNameJoin[{exportDir, "reconstruction.png"}], fig5];

(* ---------------------------------------------------------------------------
   Release gate: assert exports exist and are non-empty
   --------------------------------------------------------------------------- *)
Print["\n[Gate] verifying exports exist"];
required = {"qsp_data.json", "validation_results.json", "validation_report.md",
  "chebyshev_verification.csv", "phase_angles.csv", "eigenvalue_transform.csv",
  "chebyshev_qsp.png", "cos_qsp.png", "sin_qsp.png", "eigenvalue_transform.png",
  "reconstruction.png"};
Do[
  With[{p = FileNameJoin[{exportDir, f}]},
    If[FileExistsQ[p] && FileByteCount[p] > 0,
      Print["  OK  ", f, "  (", FileByteCount[p], " bytes)"],
      Print["  MISSING/EMPTY  ", f]]],
  {f, required}];

Print["\nBuild complete."];
