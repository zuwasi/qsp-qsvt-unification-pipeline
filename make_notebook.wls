#!/usr/bin/env wolframscript
(* ============================================================================
   make_notebook.wls  --  generate the interactive companion notebook
   QSPUnification.nb from the QSPUnification.wl package.

   Runs headlessly in a plain Wolfram kernel (no notebook front end), so it can
   never block on front-end activation or teardown. Results are computed
   in-kernel and typeset with ToBoxes (so plots embed as real GraphicsBox),
   paired with their input cells, and written directly to a proper .nb next to
   this script.

   Run:  wolframscript -file make_notebook.wls
   ============================================================================ *)

projectDir = DirectoryName[$InputFileName];
nbPath = FileNameJoin[{projectDir, "QSPUnification.nb"}];

Print["Loading package ..."];
SetDirectory[projectDir];
Get[FileNameJoin[{projectDir, "src", "QSPUnification.wl"}]];

(* ---- compute every result in-kernel (graphics -> real GraphicsBox) ---- *)
Print["Computing results ..."];
d1 = 5; phases0 = ConstantArray[0., d1 + 1];
u = QSPUnitary[phases0, 0.6];

res1 = Dataset[Table[ChebyshevDerivation[d], {d, 1, 6}]];
res2 = {Re[u[[1, 1]]], ChebyshevT[d1, 0.6]};
res3 = Plot[{QSPRealResponse[ConstantArray[0., 6], a], ChebyshevT[5, a]}, {a, -1, 1},
   PlotLegends -> {"Re P (QSP, zero phases)", "T_5(a)"},
   PlotStyle -> {Thick, Dashed}, AxesLabel -> {"a", ""},
   PlotLabel -> "Zero-phase QSP reproduces the Chebyshev polynomial"];

cosData = JacobiAngerCos[2.0, 6];
sinData = JacobiAngerSin[2.0, 6];
res4 = {"cos degree" -> cosData["Degree"], "sin degree" -> sinData["Degree"]};
res5 = Plot[cosData["Function"][x] - Cos[2.0 x], {x, -1, 1},
   PlotLabel -> "Jacobi-Anger cos truncation error (t=2, K=6)",
   AxesLabel -> {"x", "error"}];

H = BuildHamiltonian[];
res6 = Column[{
    "dimension: " <> ToString[H["Dimension"]],
    "eigenvalues: " <> ToString[H["Eigenvalues"]],
    "spectral norm: " <> ToString[Max[Abs[H["Eigenvalues"]]]]}];

eta = 0.1; t = 2.0;
cosPh = FindQSPPhases[Function[x, (1 - eta) Cos[t x]], 12, "Restarts" -> 6];
sinPh = FindQSPPhases[Function[x, (1 - eta) Sin[t x]], 13, "Restarts" -> 6];
res7 = {"cos max error" -> cosPh["MaxError"], "sin max error" -> sinPh["MaxError"]};

recon = ReconstructEvolution[H["Matrix"], t, cosPh["Phases"], sinPh["Phases"], eta];
res8 = {"max abs error" -> recon["MaxAbsError"], "overlap" -> recon["NormalizedOverlap"]};

res9 = Dataset[RunValidationSuite[]];

(* ---- cell builders ---- *)
sec[s_]  := Cell[s, "Section"];
txt[s_]  := Cell[s, "Text"];
io[in_, val_] := Sequence[
   Cell[in, "Input"],
   Cell[BoxData[ToBoxes[val]], "Output"]];

cells = {
  Cell["QSP / QSVT Unification \[LongDash] Interactive Notebook", "Title"],

  txt["Executable companion to the QSP/QSVT pipeline reproducing the core \
constructions of Martyn, Rossi, Tan & Chuang, \"Grand Unification of Quantum \
Algorithms\" (PRX Quantum 2, 040203, 2021) and the QSP Hamiltonian-simulation \
construction of Low & Chuang, \"Optimal Hamiltonian Simulation by Quantum \
Signal Processing\" (PRL 118, 010501, 2017). Outputs shown here were generated \
by make_notebook.wls; re-evaluate any cell to reproduce. Convention (Wx): \
W(a)=e^{i arccos(a) X}, S(\[Phi])=e^{i \[Phi] Z}, P(a)=\[LeftAngleBracket]0|U|0\[RightAngleBracket]."],

  sec["Setup \[LongDash] load the package"],
  txt["Loads the reusable package from this notebook's own folder."],
  Cell["SetDirectory[NotebookDirectory[]];\nGet[FileNameJoin[{NotebookDirectory[], \"src\", \"QSPUnification.wl\"}]];", "Input"],

  sec["Demo 1 \[LongDash] QSP \[DoubleLongRightArrow] Chebyshev polynomials"],
  txt["Symbolic derivation: zero-phase QSP of degree d yields P = T_d and \
Q = U_{d-1}, satisfying the Pell/Chebyshev identity T_d^2+(1-a^2)U_{d-1}^2=1."],
  io["Dataset[Table[ChebyshevDerivation[d], {d, 1, 6}]]", res1],
  txt["Numerical check at a=0.6, d=5: Re P_5(0.6) equals T_5(0.6)."],
  io["d1 = 5;\nphases0 = ConstantArray[0., d1 + 1];\nu = QSPUnitary[phases0, 0.6];\n{Re[u[[1, 1]]], ChebyshevT[d1, 0.6]}", res2],
  txt["The zero-phase QSP response traces the Chebyshev polynomial exactly."],
  io["Plot[{QSPRealResponse[ConstantArray[0., 6], a], ChebyshevT[5, a]}, {a, -1, 1},\n PlotLegends -> {\"Re P (QSP, zero phases)\", \"T_5(a)\"},\n PlotStyle -> {Thick, Dashed}, AxesLabel -> {\"a\", \"\"},\n PlotLabel -> \"Zero-phase QSP reproduces the Chebyshev polynomial\"]", res3],

  sec["Demo 2 \[LongDash] Hamiltonian simulation e^{-iHt}"],
  txt["Jacobi\[Dash]Anger expansion approximates cos(t a) and sin(t a) by \
polynomials whose coefficients are Bessel functions (t=2, K=6)."],
  io["cosData = JacobiAngerCos[2.0, 6];\nsinData = JacobiAngerSin[2.0, 6];\n{\"cos degree\" -> cosData[\"Degree\"], \"sin degree\" -> sinData[\"Degree\"]}", res4],
  io["Plot[cosData[\"Function\"][x] - Cos[2.0 x], {x, -1, 1},\n PlotLabel -> \"Jacobi-Anger cos truncation error (t=2, K=6)\",\n AxesLabel -> {\"x\", \"error\"}]", res5],
  txt["Build a driven truncated-Fock Hamiltonian from the Quantum Framework's \
second-quantization operators, normalized so the spectral norm is \[LessEqual] 1."],
  io["H = BuildHamiltonian[];\nColumn[{\n  \"dimension: \" <> ToString[H[\"Dimension\"]],\n  \"eigenvalues: \" <> ToString[H[\"Eigenvalues\"]],\n  \"spectral norm: \" <> ToString[Max[Abs[H[\"Eigenvalues\"]]]]}]", res6],
  txt["Generate symmetric QSP phase angles so Re P approximates the (scaled) \
cos and sin targets."],
  io["eta = 0.1; t = 2.0;\ncosPh = FindQSPPhases[Function[x, (1 - eta) Cos[t x]], 12, \"Restarts\" -> 6];\nsinPh = FindQSPPhases[Function[x, (1 - eta) Sin[t x]], 13, \"Restarts\" -> 6];\n{\"cos max error\" -> cosPh[\"MaxError\"], \"sin max error\" -> sinPh[\"MaxError\"]}", res7],
  txt["Reconstruct e^{-iHt} from the QSP cos/sin responses and compare to the \
exact matrix exponential."],
  io["recon = ReconstructEvolution[H[\"Matrix\"], t, cosPh[\"Phases\"], sinPh[\"Phases\"], eta];\n{\"max abs error\" -> recon[\"MaxAbsError\"], \"overlap\" -> recon[\"NormalizedOverlap\"]}", res8],

  sec["Full numerical validation suite"],
  txt["Runs the package's complete validation suite (same checks exported by \
build_project.wls)."],
  io["Dataset[RunValidationSuite[]]", res9],

  sec["References"],
  txt["[1] J. M. Martyn, Z. M. Rossi, A. K. Tan, I. L. Chuang, Grand \
Unification of Quantum Algorithms, PRX Quantum 2, 040203 (2021).\n[2] G. H. \
Low, I. L. Chuang, Optimal Hamiltonian Simulation by Quantum Signal \
Processing, Phys. Rev. Lett. 118, 010501 (2017)."]
};

Print["Writing notebook (headless kernel, no front end) ..."];
Export[nbPath, Notebook[cells, WindowTitle -> "QSP / QSVT Unification"], "NB"];

If[FileExistsQ[nbPath],
  Print["Wrote ", nbPath, "  (", FileByteCount[nbPath], " bytes)"],
  Print["ERROR: notebook not written"]; Exit[1]];
Exit[0]
