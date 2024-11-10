(* ::Package:: *)

Import["/home/string-12/sdpb_gravity/SDPB2.m"]

mu = ToExpression[ReadList["/home/string-12/sdpb_gravity/range.txt"][[1]]]
spin = 0;
bound = ToExpression[ReadList["/home/string-12/sdpb_gravity/range.txt"][[2]]];
a1 = ToExpression[ReadList["/home/string-12/sdpb_gravity/range.txt"][[3]]]

H[1,d_,J_]:=D[Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-x)/2],x]/.x-> 1

A[m_,J_,u_,d_]:=((2m^2+u)Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2u)/m^2))/2])/(m^2 (m^2+u)^2)-u^2/m^6 (((4m^2+3u)Hypergeometric2F1[-J,J+d-3,(d-2)/2,0])/(m^2+u)^2+(4 u H[1, d, J])/(m^4-u^2))
X[k_,u_,m_,J_,d_]:=(2m^2+u)/(u m^2 (m^2+u)) (m^2 Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2u)/m^2))/2])/(u m^2 (m^2+u))^(k/2)-Residue[((2m^2+t)(m^2-t)(m^2+2t))/(m^2 (u-t)t(m^2-u)(m^2+t)(m^2+u+t)) (m^2 Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2t)/m^2))/2])/(t m^2 (m^2+t))^(k/2),{t,0}]

K1[n_, J_, q_] := \[Integral]p^n A[m, J, -p^2, 7]\[DifferentialD]p   /.{p -> 1 , m -> q};
K0[n_, J_, q_] := \[Integral]p^n A[m, J, -p^2, 7]\[DifferentialD]p /.{p -> 0 , m -> q};
K[n_, J_, m_] := K1[n, J, m] - K0[n, J, m];

M1[n_, J_, q_] := \[Integral]p^n X[4,-p^2, m, J, 7]\[DifferentialD]p /.{p -> 1 , m -> q};
M0[n_, J_, q_] := \[Integral]p^n X[4,-p^2, m, J, 7]\[DifferentialD]p /.{p -> 0 , m -> q};
M[n_, J_, m_] := M1[n, J, m] - M0[n, J, m];


P1[n_, J_, q_] := \[Integral]p^n X[6,-p^2, m, J, 7]\[DifferentialD]p /.{p -> 1 , m -> q};
P0[n_, J_, q_] := \[Integral]p^n X[6,-p^2, m, J, 7]\[DifferentialD]p /.{p -> 0 , m -> q};
P[n_, J_, m_] := P1[n, J, m] - P0[n, J, m];


G[n_, D_, b_] := (2^n Gamma[(D-2)/2]Gamma[(n+1)/2])/Gamma[(D-n-3)/2] 1/b^(n+1);

B[n_, D_, b_] := (2^((D-3)/2) Gamma[(D-2)/2])/Sqrt[\[Pi]] 1/b^((D-1)/2);
F[n_, D_, b_] := (2^((D-9)/2) Gamma[(D-2)/2])/Sqrt[\[Pi]] (-27 + 12D - D^2 - 8n)/b^((D-3)/2);

(*
K[n_, J_, m_]:=NIntegrate[p^nA[m, J, -p^2, 7], {p, 0, 1}];
M[n_, J_, m_]:=NIntegrate[p^nX[4, -p^2,m, J, 7], {p, 0.02, 1}];
P[n_, J_, m_]:=NIntegrate[p^nX[6, -p^2,m, J, 7], {p, 0.02, 1}];
*)


Poly[J_, m_]:=PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, J, m], K[3, J, m], K[4, J, m], K[5, J, m], K[6, J, m], K[7, J, m], M[0, J, m], M[1, J, m], M[2, J, m], M[3, J, m], M[4, J, m], M[5, J, m], P[0, J, m], P[1, J, m], P[2, J, m], P[3, J, m]}}}];

NewPoly[b_]:=PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{G[2, 7, b] + B[2, 7, b], G[3, 7, b] + B[3, 7, b], G[4, 7, b] + B[4, 7, b], G[5, 7, b] + B[5, 7, b], G[6, 7, b] + B[6, 7, b], G[7, 7, b] + B[7, 7, b], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{F[2, 7, b], F[3, 7, b], F[4, 7, b], F[5, 7, b], F[6, 7, b], F[7, 7, b], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
        {{F[2, 7, b], F[3, 7, b], F[4, 7, b], F[5, 7, b], F[6, 7, b], F[7, 7, b], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{G[2, 7, b] - B[2, 7, b], G[3, 7, b] - B[3, 7, b], G[4, 7, b] - B[4, 7, b], G[5, 7, b] - B[5, 7, b], G[6, 7, b] - B[6, 7, b], G[7, 7, b] - B[7, 7, b], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}}];
(*        
Poly[J_, m_]:=PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, J, m], K[3, J, m], K[4, J, m], K[5, J, m], K[6, J, m], K[7, J, m]}}}];        
*)
        
TSDP[datfile_] := Module[

    {


        pols = 
        
        
       Flatten[{
       Flatten[Table[Poly[i, 1/j],{i, 0, 40, 2},{j, 1/(20*mu), 1/mu, 1/(20*mu)}]],
       {NewPoly[10], NewPoly[15], NewPoly[20], NewPoly[25], NewPoly[30], NewPoly[35], NewPoly[40]}},1],

        norm =  -1 * {19/15, 1, 29/35, 17/24, 13/21, 11/20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        obj  = -1 * ({1, 1/2, 1/3, 1/4, 1/5, 1/6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + 2 * (0)  * {1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}+0 * {1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
        (*
        norm =  {1/5, 1/6, 1/7, 1/8, 1/9, 1/10},
        obj  =  {1, 1/2, 1/3, 1/4, 1/5, 1/6}
        *)       
    },

    Print[norm];
    Print[obj];
    WritePmpJson[datfile, SDP[obj, norm, pols]]]

TSDP["out.json"]
