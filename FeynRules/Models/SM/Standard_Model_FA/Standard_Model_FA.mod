(* Patched for use with FeynCalc *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"Standard Model",
  Authors -> {"N. Christensen", "C. Duhr", "B. Fuks"},
  Version -> "1.4.7",
  Date -> "28. 09. 2016",
  Institutions -> {"Michigan State University", "Universite catholique de Louvain (CP3)", "IPHC Strasbourg / University of Strasbourg"},
  Emails -> {"neil@pa.msu.edu", "claude.duhr@uclouvain.be", "benjamin.fuks@cnrs.in2p3.fr"},
  URLs -> "http://feynrules.phys.ucl.ac.be/view/Main/StandardModel"};

FR$ClassesTranslation={};

FR$InteractionOrderPerturbativeExpansion={{QCD, 0}, {QED, 0}};

FR$GoldstoneList={S[2], S[3]};

(*     Declared indices    *)

IndexRange[ Index[Gluon] ] = NoUnfold[ Range[ 8 ] ]

IndexRange[ Index[SU2W] ] = Range[ 3 ]

IndexRange[ Index[Generation] ] = Range[ 3 ]

IndexRange[ Index[Colour] ] = NoUnfold[ Range[ 3 ] ]

IndexRange[ Index[SU2D] ] = Range[ 2 ]

(*     Declared particles    *)

M$ClassesDescription = {
V[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> "a",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> 0,
    Indices -> {} },

V[2] == {
    SelfConjugate -> True,
    PropagatorLabel -> "Z",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> FCGV["MZ"],
    Indices -> {} },

V[3] == {
    SelfConjugate -> False,
    QuantumNumbers -> {Q},
    PropagatorLabel -> "W",
    PropagatorType -> Sine,
    PropagatorArrow -> Forward,
    Mass -> FCGV["MW"],
    Indices -> {} },

V[4] == {
    SelfConjugate -> True,
    Indices -> {Index[Gluon]},
    PropagatorLabel -> "G",
    PropagatorType -> Cycles,
    PropagatorArrow -> None,
    Mass -> 0 },

U[1] == {
    SelfConjugate -> False,
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> "uA",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> 0,
    Indices -> {} },

U[2] == {
    SelfConjugate -> False,
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> "uZ",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> FCGV["MZ"],
    Indices -> {} },

U[31] == {
    SelfConjugate -> False,
    QuantumNumbers -> {GhostNumber, Q},
    PropagatorLabel -> "uWp",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> FCGV["MW"],
    Indices -> {} },

U[32] == {
    SelfConjugate -> False,
    QuantumNumbers -> {GhostNumber, -Q},
    PropagatorLabel -> "uWm",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> FCGV["MW"],
    Indices -> {} },

U[4] == {
    SelfConjugate -> False,
    Indices -> {Index[Gluon]},
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> "uG",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> 0 },

F[1] == {
    Indices -> {Index[Generation]},
    SelfConjugate -> False,
    QuantumNumbers -> {LeptonNumber},
    PropagatorLabel -> "v",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> 0 },

F[2] == {
    Indices -> {Index[Generation]},
    SelfConjugate -> False,
    QuantumNumbers -> {-Q, LeptonNumber},
    PropagatorLabel -> "l",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Ml },

F[3] == {
    Indices -> {Index[Generation], Index[Colour]},
    SelfConjugate -> False,
    QuantumNumbers -> {(2*Q)/3},
    PropagatorLabel -> "uq",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Mu },

F[4] == {
    Indices -> {Index[Generation], Index[Colour]},
    SelfConjugate -> False,
    QuantumNumbers -> {-1/3*Q},
    PropagatorLabel -> "dq",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> Md },

S[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> "H",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> FCGV["MH"],
    Indices -> {} },

S[2] == {
    SelfConjugate -> True,
    PropagatorLabel -> "Go",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> FCGV["MZ"],
    Indices -> {} },

S[3] == {
    SelfConjugate -> False,
    QuantumNumbers -> {Q},
    PropagatorLabel -> "GP",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> FCGV["MW"],
    Indices -> {} }
}


(*        Definitions       *)

FAGaugeXi[ V[1] ] = FAGaugeXi[A];
FAGaugeXi[ V[2] ] = FAGaugeXi[Z];
FAGaugeXi[ V[3] ] = FAGaugeXi[W];
FAGaugeXi[ V[4] ] = FAGaugeXi[G];
FAGaugeXi[ U[1] ] = FAGaugeXi[A];
FAGaugeXi[ U[2] ] = FAGaugeXi[Z];
FAGaugeXi[ U[31] ] = FAGaugeXi[W];
FAGaugeXi[ U[32] ] = FAGaugeXi[W];
FAGaugeXi[ U[4] ] = FAGaugeXi[G];
FAGaugeXi[ S[1] ] = 1;
FAGaugeXi[ S[2] ] = FAGaugeXi[Z];
FAGaugeXi[ S[3] ] = FAGaugeXi[W];

FCGV["MZ"][ ___ ] := FCGV["MZ"];
FCGV["MW"][ ___ ] := FCGV["MW"];
Ml[ 1 ] := Me;
Ml[ 2 ] := MMU;
Ml[ 3 ] := MTA;
Mu[ 1, _ ] := FCGV["MU"];
Mu[ 1 ] := FCGV["MU"];
Mu[ 2, _ ] := FCGV["MC"];
Mu[ 2 ] := FCGV["MC"];
Mu[ 3, _ ] := FCGV["MT"];
Mu[ 3 ] := FCGV["MT"];
Md[ 1, _ ] := FCGV["MD"];
Md[ 1 ] := FCGV["MD"];
Md[ 2, _ ] := FCGV["MS"];
Md[ 2 ] := FCGV["MS"];
Md[ 3, _ ] := FCGV["MB"];
Md[ 3 ] := FCGV["MB"];
FCGV["MH"][ ___ ] := FCGV["MH"];


TheLabel[ V[4, {__}] ] := TheLabel[V[4]];
TheLabel[ U[4, {__}] ] := TheLabel[U[4]];
TheLabel[ F[1, {1}] ] := "ve";
TheLabel[ F[1, {2}] ] := "vm";
TheLabel[ F[1, {3}] ] := "vt";
TheLabel[ F[2, {1}] ] := "e";
TheLabel[ F[2, {2}] ] := "mu";
TheLabel[ F[2, {3}] ] := "ta";
TheLabel[ F[3, {1, _}] ] := "u";
TheLabel[ F[3, {1}] ] := "u";
TheLabel[ F[3, {2, _}] ] := "c";
TheLabel[ F[3, {2}] ] := "c";
TheLabel[ F[3, {3, _}] ] := "t";
TheLabel[ F[3, {3}] ] := "t";
TheLabel[ F[4, {1, _}] ] := "d";
TheLabel[ F[4, {1}] ] := "d";
TheLabel[ F[4, {2, _}] ] := "s";
TheLabel[ F[4, {2}] ] := "s";
TheLabel[ F[4, {3, _}] ] := "b";
TheLabel[ F[4, {3}] ] := "b";


(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[2] , S[2] , S[2] , S[2] ] == {{(-6*I)*lam, 0}},

C[ S[2] , S[2] , S[3] , -S[3] ] == {{(-2*I)*lam, 0}},

C[ S[3] , S[3] , -S[3] , -S[3] ] == {{(-4*I)*lam, 0}},

C[ S[2] , S[2] , S[1] , S[1] ] == {{(-2*I)*lam, 0}},

C[ S[3] , -S[3] , S[1] , S[1] ] == {{(-2*I)*lam, 0}},

C[ S[1] , S[1] , S[1] , S[1] ] == {{(-6*I)*lam, 0}},

C[ S[2] , S[2] , S[1] ] == {{(-2*I)*lam*vev, 0}},

C[ S[3] , -S[3] , S[1] ] == {{(-2*I)*lam*vev, 0}},

C[ S[1] , S[1] , S[1] ] == {{(-6*I)*lam*vev, 0}},

C[ S[3] , -S[3] , V[1] , V[1] ] == {{(2*I)*FCGV["EL"]^2, 0}},

C[ S[3] , -S[3] , V[1] ] == {{(-I)*gc11, 0}, {I*gc11, 0}},

C[ -U[1] , U[32] , V[3] ] == {{I*gc12, 0}, {I*gc12, 0}, {0, 0}},

C[ -U[1] , U[31] , -V[3] ] == {{I*gc13, 0}, {I*gc13, 0}, {0, 0}},

C[ -S[3] , -U[32] , U[1] ] == {{(FCGV["EL"]^2*vev)/(2*sw), 0}},

C[ -U[32] , U[1] , -V[3] ] == {{I*gc15, 0}, {I*gc15, 0}, {0, 0}},

C[ S[2] , -U[32] , U[32] ] == {{-1/4*(FCGV["EL"]^2*vev)/sw^2, 0}},

C[ S[1] , -U[32] , U[32] ] == {{((-1/4*I)*FCGV["EL"]^2*vev)/sw^2, 0}},

C[ -U[32] , U[32] , V[1] ] == {{I*gc18, 0}, {I*gc18, 0}, {0, 0}},

C[ -U[32] , U[32] , V[2] ] == {{I*gc19, 0}, {I*gc19, 0}, {0, 0}},

C[ -S[3] , -U[32] , U[2] ] == {{(FCGV["EL"]^2*(cw - sw)*(cw + sw)*vev)/(4*cw*sw^2), 0}},

C[ -U[32] , U[2] , -V[3] ] == {{I*gc21, 0}, {I*gc21, 0}, {0, 0}},

C[ S[3] , -U[31] , U[1] ] == {{-1/2*(FCGV["EL"]^2*vev)/sw, 0}},

C[ -U[31] , U[1] , V[3] ] == {{I*gc23, 0}, {I*gc23, 0}, {0, 0}},

C[ S[2] , -U[31] , U[31] ] == {{(FCGV["EL"]^2*vev)/(4*sw^2), 0}},

C[ S[1] , -U[31] , U[31] ] == {{((-1/4*I)*FCGV["EL"]^2*vev)/sw^2, 0}},

C[ -U[31] , U[31] , V[1] ] == {{I*gc26, 0}, {I*gc26, 0}, {0, 0}},

C[ -U[31] , U[31] , V[2] ] == {{I*gc27, 0}, {I*gc27, 0}, {0, 0}},

C[ S[3] , -U[31] , U[2] ] == {{-1/4*(FCGV["EL"]^2*(cw - sw)*(cw + sw)*vev)/(cw*sw^2), 0}},

C[ -U[31] , U[2] , V[3] ] == {{I*gc29, 0}, {I*gc29, 0}, {0, 0}},

C[ S[3] , -U[2] , U[32] ] == {{(FCGV["EL"]^2*(cw^2 + sw^2)*vev)/(4*cw*sw^2), 0}},

C[ -U[2] , U[32] , V[3] ] == {{I*gc31, 0}, {I*gc31, 0}, {0, 0}},

C[ -S[3] , -U[2] , U[31] ] == {{-1/4*(FCGV["EL"]^2*(cw^2 + sw^2)*vev)/(cw*sw^2), 0}},

C[ -U[2] , U[31] , -V[3] ] == {{I*gc33, 0}, {I*gc33, 0}, {0, 0}},

C[ S[1] , -U[2] , U[2] ] == {{((-1/4*I)*FCGV["EL"]^2*(cw^2 + sw^2)^2*vev)/(cw^2*sw^2), 0}},

C[ -U[4, {e1x1}] , U[4, {e2x1}] , V[4, {e3x2}] ] == {{gc35*FASUNF[e3x2, e1x1, e2x1], 0}, {gc35*FASUNF[e3x2, e1x1, e2x1], 0}, {0, 0}},

C[ V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] ] == {{-(gc36*FASUNF[e1x2, e2x2, e3x2]), 0}, {gc36*FASUNF[e1x2, e2x2, e3x2], 0}, {gc36*FASUNF[e1x2, e2x2, e3x2], 0}, {-(gc36*FASUNF[e1x2, e2x2, e3x2]), 0}, {-(gc36*FASUNF[e1x2, e2x2, e3x2]), 0}, {gc36*FASUNF[e1x2, e2x2, e3x2], 0}},

C[ V[4, {e1x2}] , V[4, {e2x2}] , V[4, {e3x2}] , V[4, {e4x2}] ] == {{(-I)*gc37*(FASUNF[e1x2, e2x2, e3x2, e4x2] + FASUNF[e1x2, e3x2, e2x2, e4x2]), 0}, {I*gc37*(FASUNF[e1x2, e2x2, e3x2, e4x2] - FASUNF[e1x2, e4x2, e2x2, e3x2]), 0}, {I*gc37*(FASUNF[e1x2, e3x2, e2x2, e4x2] + FASUNF[e1x2, e4x2, e2x2, e3x2]), 0}},

C[ -F[4, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , -S[3] ] == {{gc38L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {gc38R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , S[2] ] == {{gc39L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {gc39R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , S[1] ] == {{I*gc40L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc40R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[2, {e1x2}] , F[1, {e2x2}] , -S[3] ] == {{gc41[e1x2, e2x2], 0}, {0, 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , S[2] ] == {{gc42L[e1x2, e2x2], 0}, {gc42R[e1x2, e2x2], 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , S[1] ] == {{I*gc43L[e1x2, e2x2], 0}, {I*gc43R[e1x2, e2x2], 0}},

C[ -F[3, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , S[3] ] == {{gc44L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {gc44R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , S[2] ] == {{gc45L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {gc45R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , S[1] ] == {{I*gc46L[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc46R[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ S[2] , -S[3] , V[1] , V[3] ] == {{((-1/2*I)*FCGV["EL"]^2)/sw, 0}},

C[ -S[3] , S[1] , V[1] , V[3] ] == {{-1/2*FCGV["EL"]^2/sw, 0}},

C[ -S[3] , V[1] , V[3] ] == {{-1/2*(FCGV["EL"]^2*vev)/sw, 0}},

C[ S[2] , -S[3] , V[3] ] == {{(-I)*gc50, 0}, {I*gc50, 0}},

C[ -S[3] , S[1] , V[3] ] == {{-gc51, 0}, {gc51, 0}},

C[ V[1] , V[3] , -V[3] ] == {{(-I)*gc52, 0}, {I*gc52, 0}, {I*gc52, 0}, {(-I)*gc52, 0}, {(-I)*gc52, 0}, {I*gc52, 0}},

C[ S[2] , S[3] , V[1] , -V[3] ] == {{((-1/2*I)*FCGV["EL"]^2)/sw, 0}},

C[ S[3] , S[1] , V[1] , -V[3] ] == {{FCGV["EL"]^2/(2*sw), 0}},

C[ S[3] , V[1] , -V[3] ] == {{(FCGV["EL"]^2*vev)/(2*sw), 0}},

C[ S[2] , S[3] , -V[3] ] == {{(-I)*gc56, 0}, {I*gc56, 0}},

C[ S[3] , S[1] , -V[3] ] == {{-gc57, 0}, {gc57, 0}},

C[ S[2] , S[2] , V[3] , -V[3] ] == {{((I/2)*FCGV["EL"]^2)/sw^2, 0}},

C[ S[3] , -S[3] , V[3] , -V[3] ] == {{((I/2)*FCGV["EL"]^2)/sw^2, 0}},

C[ S[1] , S[1] , V[3] , -V[3] ] == {{((I/2)*FCGV["EL"]^2)/sw^2, 0}},

C[ S[1] , V[3] , -V[3] ] == {{((I/2)*FCGV["EL"]^2*vev)/sw^2, 0}},

C[ V[1] , V[1] , V[3] , -V[3] ] == {{(-I)*gc62, 0}, {(-I)*gc62, 0}, {(2*I)*gc62, 0}},

C[ V[3] , -V[3] , V[2] ] == {{(-I)*gc63, 0}, {I*gc63, 0}, {I*gc63, 0}, {(-I)*gc63, 0}, {(-I)*gc63, 0}, {I*gc63, 0}},

C[ V[3] , V[3] , -V[3] , -V[3] ] == {{(-I)*gc64, 0}, {(-I)*gc64, 0}, {(2*I)*gc64, 0}},

C[ -F[1, {e1x2}] , F[2, {e2x2}] , S[3] ] == {{0, 0}, {gc65R[e1x2, e2x2], 0}},

C[ S[3] , -S[3] , V[1] , V[2] ] == {{(I*FCGV["EL"]^2*(cw - sw)*(cw + sw))/(cw*sw), 0}},

C[ S[2] , S[1] , V[2] ] == {{-gc67, 0}, {gc67, 0}},

C[ S[3] , -S[3] , V[2] ] == {{(-I)*gc68, 0}, {I*gc68, 0}},

C[ S[2] , -S[3] , V[3] , V[2] ] == {{((I/2)*FCGV["EL"]^2)/cw, 0}},

C[ -S[3] , S[1] , V[3] , V[2] ] == {{FCGV["EL"]^2/(2*cw), 0}},

C[ -S[3] , V[3] , V[2] ] == {{(FCGV["EL"]^2*vev)/(2*cw), 0}},

C[ S[2] , S[3] , -V[3] , V[2] ] == {{((I/2)*FCGV["EL"]^2)/cw, 0}},

C[ S[3] , S[1] , -V[3] , V[2] ] == {{-1/2*FCGV["EL"]^2/cw, 0}},

C[ S[3] , -V[3] , V[2] ] == {{-1/2*(FCGV["EL"]^2*vev)/cw, 0}},

C[ V[1] , V[3] , -V[3] , V[2] ] == {{(-2*I)*gc75, 0}, {I*gc75, 0}, {I*gc75, 0}},

C[ S[2] , S[2] , V[2] , V[2] ] == {{((I/2)*FCGV["EL"]^2*(cw^2 + sw^2)^2)/(cw^2*sw^2), 0}},

C[ S[3] , -S[3] , V[2] , V[2] ] == {{((I/2)*FCGV["EL"]^2*(cw - sw)^2*(cw + sw)^2)/(cw^2*sw^2), 0}},

C[ S[1] , S[1] , V[2] , V[2] ] == {{((I/2)*FCGV["EL"]^2*(cw^2 + sw^2)^2)/(cw^2*sw^2), 0}},

C[ S[1] , V[2] , V[2] ] == {{((I/2)*FCGV["EL"]^2*(cw^2 + sw^2)^2*vev)/(cw^2*sw^2), 0}},

C[ V[3] , -V[3] , V[2] , V[2] ] == {{(-I)*gc80, 0}, {(-I)*gc80, 0}, {(2*I)*gc80, 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , V[1] ] == {{I*gc81*IndexDelta[e1x2, e2x2], 0}, {I*gc81*IndexDelta[e1x2, e2x2], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[1] ] == {{I*gc82*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc82*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[1] ] == {{I*gc83*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc83*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[4, {e3x2}] ] == {{I*gc84*IndexDelta[e1x2, e2x2]*FASUNT[e3x2, e1x3, e2x3], 0}, {I*gc84*IndexDelta[e1x2, e2x2]*FASUNT[e3x2, e1x3, e2x3], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[4, {e3x2}] ] == {{I*gc85*IndexDelta[e1x2, e2x2]*FASUNT[e3x2, e1x3, e2x3], 0}, {I*gc85*IndexDelta[e1x2, e2x2]*FASUNT[e3x2, e1x3, e2x3], 0}},

C[ -F[3, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[3] ] == {{I*gc86[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {0, 0}},

C[ -F[4, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , -V[3] ] == {{I*gc87[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {0, 0}},

C[ -F[1, {e1x2}] , F[2, {e2x2}] , V[3] ] == {{I*gc88*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -F[2, {e1x2}] , F[1, {e2x2}] , -V[3] ] == {{I*gc89*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -F[3, {e1x2, e1x3}] , F[3, {e2x2, e2x3}] , V[2] ] == {{I*gc90L*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc90R*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[4, {e1x2, e1x3}] , F[4, {e2x2, e2x3}] , V[2] ] == {{I*gc91L*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}, {I*gc91R*IndexDelta[e1x2, e2x2]*IndexDelta[e1x3, e2x3], 0}},

C[ -F[1, {e1x2}] , F[1, {e2x2}] , V[2] ] == {{I*gc92*IndexDelta[e1x2, e2x2], 0}, {0, 0}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , V[2] ] == {{I*gc93L*IndexDelta[e1x2, e2x2], 0}, {I*gc93R*IndexDelta[e1x2, e2x2], 0}}

}

(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)

(* Parameter replacement lists (These lists were created by FeynRules) *)

(* FA Couplings *)

M$FACouplings = {
     gc11 -> -FCGV["EL"],
     gc12 -> FCGV["EL"],
     gc13 -> -FCGV["EL"],
     gc15 -> FCGV["EL"],
     gc18 -> -FCGV["EL"],
     gc19 -> -((cw*FCGV["EL"])/sw),
     gc21 -> (cw*FCGV["EL"])/sw,
     gc23 -> -FCGV["EL"],
     gc26 -> FCGV["EL"],
     gc27 -> (cw*FCGV["EL"])/sw,
     gc29 -> -((cw*FCGV["EL"])/sw),
     gc31 -> (cw*FCGV["EL"])/sw,
     gc33 -> -((cw*FCGV["EL"])/sw),
     gc35 -> FAGS,
     gc36 -> -FAGS,
     gc37 -> -FAGS^2,
     gc38L[e1x2_, e2x2_] -> IndexSum[Conjugate[CKM[e2x2, Generation$1]]*Conjugate[yd[Generation$1, e1x2]], {Generation$1, 1, 3}],
     gc38R[e1x2_, e2x2_] -> -IndexSum[Conjugate[CKM[Generation$1, e1x2]]*yu[Generation$1, e2x2], {Generation$1, 1, 3}],
     gc39L[e1x2_, e2x2_] -> -(Conjugate[yd[e2x2, e1x2]]/Sqrt[2]),
     gc39R[e1x2_, e2x2_] -> yd[e1x2, e2x2]/Sqrt[2],
     gc40L[e1x2_, e2x2_] -> -(Conjugate[yd[e2x2, e1x2]]/Sqrt[2]),
     gc40R[e1x2_, e2x2_] -> -(yd[e1x2, e2x2]/Sqrt[2]),
     gc41[e1x2_, e2x2_] -> Conjugate[yl[e2x2, e1x2]],
     gc42L[e1x2_, e2x2_] -> -(Conjugate[yl[e2x2, e1x2]]/Sqrt[2]),
     gc42R[e1x2_, e2x2_] -> yl[e1x2, e2x2]/Sqrt[2],
     gc43L[e1x2_, e2x2_] -> -(Conjugate[yl[e2x2, e1x2]]/Sqrt[2]),
     gc43R[e1x2_, e2x2_] -> -(yl[e1x2, e2x2]/Sqrt[2]),
     gc44L[e1x2_, e2x2_] -> IndexSum[CKM[Generation$1, e2x2]*Conjugate[yu[Generation$1, e1x2]], {Generation$1, 1, 3}],
     gc44R[e1x2_, e2x2_] -> -IndexSum[CKM[e1x2, Generation$1]*yd[Generation$1, e2x2], {Generation$1, 1, 3}],
     gc45L[e1x2_, e2x2_] -> Conjugate[yu[e2x2, e1x2]]/Sqrt[2],
     gc45R[e1x2_, e2x2_] -> -(yu[e1x2, e2x2]/Sqrt[2]),
     gc46L[e1x2_, e2x2_] -> -(Conjugate[yu[e2x2, e1x2]]/Sqrt[2]),
     gc46R[e1x2_, e2x2_] -> -(yu[e1x2, e2x2]/Sqrt[2]),
     gc50 -> FCGV["EL"]/(2*sw),
     gc51 -> -1/2*FCGV["EL"]/sw,
     gc52 -> FCGV["EL"],
     gc56 -> -1/2*FCGV["EL"]/sw,
     gc57 -> -1/2*FCGV["EL"]/sw,
     gc62 -> -FCGV["EL"]^2,
     gc63 -> (cw*FCGV["EL"])/sw,
     gc64 -> FCGV["EL"]^2/sw^2,
     gc65R[e1x2_, e2x2_] -> -yl[e1x2, e2x2],
     gc67 -> -1/2*(FCGV["EL"]*(cw^2 + sw^2))/(cw*sw),
     gc68 -> -1/2*(cw*FCGV["EL"])/sw + (FCGV["EL"]*sw)/(2*cw),
     gc75 -> (cw*FCGV["EL"]^2)/sw,
     gc80 -> -((cw^2*FCGV["EL"]^2)/sw^2),
     gc81 -> -FCGV["EL"],
     gc82 -> (2*FCGV["EL"])/3,
     gc83 -> -1/3*FCGV["EL"],
     gc84 -> FAGS,
     gc85 -> FAGS,
     gc86[e1x2_, e2x2_] -> (FCGV["EL"]*CKM[e1x2, e2x2])/(Sqrt[2]*sw),
     gc87[e1x2_, e2x2_] -> (FCGV["EL"]*Conjugate[CKM[e2x2, e1x2]])/(Sqrt[2]*sw),
     gc88 -> FCGV["EL"]/(Sqrt[2]*sw),
     gc89 -> FCGV["EL"]/(Sqrt[2]*sw),
     gc90L -> (cw*FCGV["EL"])/(2*sw) - (FCGV["EL"]*sw)/(6*cw),
     gc90R -> (-2*FCGV["EL"]*sw)/(3*cw),
     gc91L -> -1/6*(FCGV["EL"]*(3*cw^2 + sw^2))/(cw*sw),
     gc91R -> (FCGV["EL"]*sw)/(3*cw),
     gc92 -> (FCGV["EL"]*(cw^2 + sw^2))/(2*cw*sw),
     gc93L -> -1/2*(FCGV["EL"]*(cw^2 - sw^2))/(cw*sw),
     gc93R -> (FCGV["EL"]*sw)/cw};

