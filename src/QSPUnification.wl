(* ::Package:: *)

(* QSPUnification.wl
   Reusable Wolfram Language package reproducing the core constructions of
   "Grand Unification of Quantum Algorithms" (Martyn, Rossi, Tan, Chuang, 2021)
   and "Optimal Hamiltonian Simulation by Quantum Signal Processing"
   (Low & Chuang, 2017).

   Convention (paper's Wx convention):
     Signal operator   W(a) = e^{i theta X},  a = Cos[theta] in [-1,1]
                            = {{a, I Sqrt[1-a^2]},{I Sqrt[1-a^2], a}}
     Phase operator    S(phi) = e^{i phi Z} = {{e^{I phi},0},{0,e^{-I phi}}}
     QSP sequence      U_Phi(a) = S(phi0) . Prod_{k=1}^d ( W(a) . S(phi_k) )
                       (d signal calls, d+1 phases phi0..phi_d)
     QSP polynomial    P(a) = (U_Phi(a))_{1,1}
*)

BeginPackage["QSPUnification`"];

QSPSignalOperator::usage    = "QSPSignalOperator[a] returns the 2x2 QSP signal operator W(a)=e^{i arccos(a) X}.";
QSPPhaseOperator::usage     = "QSPPhaseOperator[phi] returns the 2x2 QSP phase operator S(phi)=e^{i phi Z}.";
QSPUnitary::usage           = "QSPUnitary[phases, a] returns the 2x2 QSP unitary U_Phi(a) for a phase list {phi0,...,phid}.";
QSPPolynomial::usage        = "QSPPolynomial[phases, a] returns the QSP polynomial P(a) = (U_Phi(a))_{1,1}.";
QSPRealResponse::usage      = "QSPRealResponse[phases, a] returns Re[P(a)], the achievable real polynomial response.";
ChebyshevDerivation::usage  = "ChebyshevDerivation[d] symbolically derives that zero-phase QSP of degree d yields T_d and U_{d-1}, and checks the Pell/Chebyshev identity. Returns an Association.";
JacobiAngerCos::usage       = "JacobiAngerCos[t, K] returns <|\"Function\"->pure fn, \"Coefficients\"->cheb coeffs, \"Degree\"->2K|> approximating Cos[t x] via the Jacobi-Anger expansion truncated at 2K.";
JacobiAngerSin::usage       = "JacobiAngerSin[t, K] returns the analogous data approximating Sin[t x] truncated at degree 2K+1.";
FindQSPPhases::usage        = "FindQSPPhases[f, d] finds symmetric QSP phases {phi0..phid} so Re[P(a)] approximates f[a] on [-1,1]. Returns <|\"Phases\"->..,\"Residual\"->..,\"MaxError\"->..|>.";
BuildHamiltonian::usage     = "BuildHamiltonian[] builds a two-mode Hamiltonian from second-quantization (creation/annihilation) operators of the Wolfram Quantum Framework, normalized so ||H||<=1, and returns an Association with matrix and spectral data.";
ReconstructEvolution::usage = "ReconstructEvolution[H, t, cosData, sinData, eta] reconstructs e^{-I H t} from QSP cos/sin responses and returns error metrics.";
RunValidationSuite::usage   = "RunValidationSuite[] runs the full numerical validation suite and returns an Association of results.";
ComplexToPairs::usage       = "ComplexToPairs[expr] serializes complex numbers/arrays to {re,im} pairs for JSON export.";

Begin["`Private`"];

(* ---------- Core QSP objects ---------- *)

pauliX = {{0, 1}, {1, 0}};
pauliZ = {{1, 0}, {0, -1}};

QSPSignalOperator[a_] := With[{th = ArcCos[a]},
  {{Cos[th], I Sin[th]}, {I Sin[th], Cos[th]}}
];

QSPPhaseOperator[phi_] := {{Exp[I phi], 0}, {0, Exp[-I phi]}};

QSPUnitary[phases_List, a_] := Module[{d, w, u},
  d = Length[phases] - 1;
  w = QSPSignalOperator[a];
  u = QSPPhaseOperator[phases[[1]]];
  Do[u = u . w . QSPPhaseOperator[phases[[k + 1]]], {k, 1, d}];
  u
];

QSPPolynomial[phases_List, a_] := QSPUnitary[phases, a][[1, 1]];
QSPRealResponse[phases_List, a_] := Re[QSPPolynomial[phases, a]];

(* ---------- Symbolic Chebyshev derivation ---------- *)

ChebyshevDerivation[d_Integer] := Module[
  {th, wth, wpow, p11, p12, tCheck, uCheck, idCheck, a},
  wth = {{Cos[th], I Sin[th]}, {I Sin[th], Cos[th]}};
  wpow = MatrixPower[wth, d];
  p11 = FullSimplify[wpow[[1, 1]]];
  p12 = FullSimplify[wpow[[1, 2]]];
  (* (1,1) entry should be Cos[d th] = T_d(Cos th) *)
  tCheck = FullSimplify[p11 - ChebyshevT[d, Cos[th]]] === 0;
  (* (1,2) entry should be I Sin[d th] = I Sin[th] U_{d-1}(Cos th) *)
  uCheck = FullSimplify[p12 - I Sin[th] ChebyshevU[d - 1, Cos[th]]] === 0;
  (* Chebyshev "Pell" identity T_d^2 + (1-a^2) U_{d-1}^2 = 1 *)
  idCheck = FullSimplify[ChebyshevT[d, a]^2 + (1 - a^2) ChebyshevU[d - 1, a]^2 - 1] === 0;
  <|
    "Degree" -> d,
    "P11" -> p11,
    "P12" -> p12,
    "P_equals_Td" -> tCheck,
    "Q_equals_Udm1" -> uCheck,
    "ChebyshevIdentity" -> idCheck
  |>
];

(* ---------- Jacobi-Anger expansions ---------- *)

JacobiAngerCos[t_, K_Integer] := Module[{coeffs, fn, x},
  (* Cos[t x] = J0(t) + 2 Sum_{k>=1} (-1)^k J_{2k}(t) T_{2k}(x) *)
  coeffs = Prepend[
    Table[2 (-1)^k BesselJ[2 k, t], {k, 1, K}],
    BesselJ[0, t]
  ];  (* coeffs for T_0, T_2, T_4, ..., T_{2K} *)
  fn = Function[x, Evaluate[
    BesselJ[0, t] + Sum[2 (-1)^k BesselJ[2 k, t] ChebyshevT[2 k, x], {k, 1, K}]
  ]];
  <|"Function" -> fn, "Coefficients" -> N[coeffs], "Orders" -> Table[2 k, {k, 0, K}],
    "Degree" -> 2 K, "Parity" -> "Even"|>
];

JacobiAngerSin[t_, K_Integer] := Module[{coeffs, fn, x},
  (* Sin[t x] = 2 Sum_{k>=0} (-1)^k J_{2k+1}(t) T_{2k+1}(x) *)
  coeffs = Table[2 (-1)^k BesselJ[2 k + 1, t], {k, 0, K}];
  fn = Function[x, Evaluate[
    Sum[2 (-1)^k BesselJ[2 k + 1, t] ChebyshevT[2 k + 1, x], {k, 0, K}]
  ]];
  <|"Function" -> fn, "Coefficients" -> N[coeffs], "Orders" -> Table[2 k + 1, {k, 0, K}],
    "Degree" -> 2 K + 1, "Parity" -> "Odd"|>
];

(* ---------- Phase-angle generation (optimization, symmetric phases) ---------- *)

Options[FindQSPPhases] = {"Nodes" -> Automatic, "Restarts" -> 6, "Tolerance" -> 1.*^-8};

FindQSPPhases[f_, d_Integer, OptionsPattern[]] := Module[
  {nPhase, nFree, M, nodes, nodeVals, target, uNum, obj, vars, start,
   best, bestVal, sol, phasesFull, resid, maxErr, seeds, restarts, tol},
  nPhase = d + 1;
  nFree = Ceiling[nPhase/2];               (* symmetric: phi_k = phi_{d-k} *)
  restarts = OptionValue["Restarts"];
  tol = OptionValue["Tolerance"];
  M = OptionValue["Nodes"] /. Automatic -> Max[2 (d + 2), 40];
  nodes = N@Table[Cos[(2 j - 1) Pi/(2 M)], {j, 1, M}];   (* Chebyshev nodes *)
  target = f /@ nodes;

  (* Build symmetric phase list from free variables (numeric-guarded) *)
  symPhases[xv_?(VectorQ[#, NumericQ] &)] :=
    Table[xv[[Min[k, nPhase + 1 - k]]], {k, 1, nPhase}];
  (* Numeric-only least-squares objective; guard keeps FindMinimum numeric *)
  objFun[xv_?(VectorQ[#, NumericQ] &)] :=
    Total[((QSPRealResponse[symPhases[xv], #] & /@ nodes) - target)^2];

  best = None; bestVal = Infinity;
  vars = Array[xx, nFree];
  SeedRandom[1234];                       (* reproducible restarts *)
  seeds = Join[
    {ConstantArray[0., nFree]},
    Table[RandomReal[{-0.4, 0.4}, nFree], {restarts - 1}]
  ];

  (* Optimization loop over restarts, keep the best minimizer.
     NB: do NOT trap FindMinimum via Check[] -- benign warnings such as
     FindMinimum::lstol would be misread as failures. Validate the shape. *)
  Do[
    Module[{res, val},
      res = Quiet@FindMinimum[objFun[vars], Evaluate[Transpose[{vars, seed}]],
          MaxIterations -> 800, Method -> "QuasiNewton",
          AccuracyGoal -> 10, PrecisionGoal -> 10];
      If[MatchQ[res, {_?NumericQ, {__Rule}}],
        val = res[[1]];
        If[val < bestVal, bestVal = val; best = vars /. res[[2]]];
      ];
    ],
    {seed, seeds}
  ];

  If[best === None, Return[<|"Phases" -> $Failed, "Residual" -> Infinity|>]];

  (* Wrap phases into (-pi, pi]; e^{i phi Z} has period 2 pi so U is unchanged *)
  best = Mod[best + Pi, 2 Pi] - Pi;
  phasesFull = symPhases[best];
  resid = Sqrt[bestVal/M];
  maxErr = Max@Abs@Table[QSPRealResponse[phasesFull, nodes[[j]]] - target[[j]], {j, 1, M}];
  <|
    "Phases" -> phasesFull,
    "FreePhases" -> best,
    "Degree" -> d,
    "Nodes" -> M,
    "Residual" -> resid,
    "MaxError" -> maxErr
  |>
];

(* ---------- Hamiltonian from second quantization ---------- *)

Options[BuildHamiltonian] = {"FockSize" -> 4, "Omega" -> 1.0, "Drive" -> 0.6, "NormTarget" -> 0.85};

BuildHamiltonian[OptionsPattern[]] := Module[
  {annFn, size, omega, g, normTarget, aOp, aMat, adagMat, nMat, mat,
   shift, Hcentered, scale, vals, vecs, Hn},
  Needs["Wolfram`QuantumFramework`"];
  size = OptionValue["FockSize"]; omega = OptionValue["Omega"];
  g = OptionValue["Drive"]; normTarget = OptionValue["NormTarget"];
  (* Genuine second-quantization ladder operator from the Wolfram Quantum
     Framework's SecondQuantization module (truncated Fock space). *)
  annFn = PacletSymbol["Wolfram/QuantumFramework", "AnnihilationOperator"];
  aOp = annFn[size, {1}];
  aMat = Normal[aOp["Matrix"]];        (* Band sqrt(1..size-1): a in Fock basis *)
  adagMat = ConjugateTranspose[aMat];  (* creation a^dagger *)
  nMat = adagMat . aMat;               (* number operator n = a^dagger a = diag(0..size-1) *)
  (* Driven harmonic oscillator:  H_phys = omega n + g (a + a^dagger) *)
  mat = omega nMat + g (aMat + adagMat);
  mat = (mat + ConjugateTranspose[mat])/2;           (* enforce Hermiticity *)
  shift = Tr[mat]/Length[mat];                        (* center spectrum (traceless) *)
  Hcentered = mat - shift IdentityMatrix[Length[mat]];
  scale = Max[Abs[Eigenvalues[Hcentered]]]/normTarget;   (* ||H|| = normTarget < 1 *)
  Hn = Hcentered/scale;
  {vals, vecs} = Eigensystem[Hn];
  With[{ord = Ordering[Re[vals]]},        (* sort eigenvalues & vectors together *)
    vals = vals[[ord]]; vecs = vecs[[ord]]];
  <|
    "Matrix" -> Hn,
    "RawMatrix" -> mat,
    "AnnihilationMatrix" -> aMat,
    "NumberMatrix" -> nMat,
    "Shift" -> shift,
    "Scale" -> scale,
    "Params" -> <|"FockSize" -> size, "Omega" -> omega, "Drive" -> g, "NormTarget" -> normTarget|>,
    "Eigenvalues" -> Re[vals],            (* sorted, index-matched to Eigenvectors *)
    "Eigenvectors" -> vecs,               (* row j is the eigenvector for Eigenvalues[[j]] *)
    "Dimension" -> Length[Hn]
  |>
];

(* ---------- Evolution reconstruction from QSP responses ---------- *)

ReconstructEvolution[H_, t_, cosPhases_, sinPhases_, eta_] := Module[
  {vals, vecs, dim, exact, recon, j, lam, pc, ps, proj, err, fid},
  dim = Length[H];
  {vals, vecs} = Eigensystem[H];
  exact = MatrixExp[-I H t];
  recon = ConstantArray[0. + 0. I, {dim, dim}];
  Do[
    lam = Re[vals[[j]]];
    pc = QSPRealResponse[cosPhases, lam]/(1 - eta);   (* approx Cos[lam t] *)
    ps = QSPRealResponse[sinPhases, lam]/(1 - eta);   (* approx Sin[lam t] *)
    proj = Outer[Times, vecs[[j]], Conjugate[vecs[[j]]]];
    recon += (pc - I ps) proj;
    , {j, 1, dim}];
  err = Max@Abs@Flatten[recon - exact];
  fid = Abs[Tr[ConjugateTranspose[exact] . recon]]/dim;
  <|
    "Exact" -> exact,
    "Reconstructed" -> recon,
    "MaxAbsError" -> err,
    "NormalizedOverlap" -> fid
  |>
];

(* ---------- Serialization helper ---------- *)

ComplexToPairs[z_?NumericQ] := {Re[N[z]], Im[N[z]]};
ComplexToPairs[arr_List] := Map[ComplexToPairs, arr];

(* ---------- Validation suite ---------- *)

RunValidationSuite[] := Module[
  {results, aTest, d, ja, jaErr, t, eta, K, h, cosD, sinD, recon},
  results = <||>;
  SeedRandom[7];
  aTest = RandomReal[{-1, 1}, 25];

  (* 1. Symbolic Chebyshev derivation for d=1..6 *)
  results["ChebyshevDerivation"] = Association@Table[
    d -> ChebyshevDerivation[d], {d, 1, 6}];

  (* 2. Numeric: zero-phase QSP Re(0,0) == T_d(a) *)
  results["ZeroPhaseChebyshev"] = Association@Table[
    d -> Max@Abs@Table[
        QSPRealResponse[ConstantArray[0., d + 1], a] - ChebyshevT[d, a], {a, aTest}],
    {d, 1, 6}];

  (* 3. QSP unitarity / Pell identity |P|^2 + (1-a^2)|Q|^2 == 1 (zero-phase, d=5) *)
  results["QSPUnitarity"] = Max@Abs@Table[
    Module[{u = QSPUnitary[ConstantArray[0., 6], a]},
      Abs[u[[1, 1]]]^2 + Abs[u[[1, 2]]]^2 - 1], {a, aTest}];

  (* 4. Jacobi-Anger truncation error vs exact cos/sin *)
  t = 2.0; K = 6;
  ja = JacobiAngerCos[t, K];
  results["JacobiAngerCosError"] = Max@Table[Abs[ja["Function"][a] - Cos[t a]], {a, aTest}];
  ja = JacobiAngerSin[t, K];
  results["JacobiAngerSinError"] = Max@Table[Abs[ja["Function"][a] - Sin[t a]], {a, aTest}];

  (* 5. Phase-angle generation accuracy for the cos/sin targets *)
  eta = 0.1;
  cosD = FindQSPPhases[Function[x, (1 - eta) Cos[t x]], 2 K, "Restarts" -> 8];
  sinD = FindQSPPhases[Function[x, (1 - eta) Sin[t x]], 2 K + 1, "Restarts" -> 8];
  results["PhaseGenCosMaxError"] = cosD["MaxError"];
  results["PhaseGenSinMaxError"] = sinD["MaxError"];

  (* 6. Full Hamiltonian-simulation reconstruction error *)
  h = BuildHamiltonian[];
  recon = ReconstructEvolution[h["Matrix"], t, cosD["Phases"], sinD["Phases"], eta];
  results["HamiltonianSimError"] = recon["MaxAbsError"];
  results["HamiltonianSimOverlap"] = recon["NormalizedOverlap"];

  results
];

End[];
EndPackage[];
