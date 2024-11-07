(* ::Package:: *)

<<"SDPB.m";

H[1,d_,J_]:=D[Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-x)/2],x]/.x-> 1

A[m_,J_,u_,d_]:=((2m^2+u)Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2u)/m^2))/2])/(m^2 (m^2+u)^2)-u^2/m^6 (((4m^2+3u)Hypergeometric2F1[-J,J+d-3,(d-2)/2,0])/(m^2+u)^2+(4 u H[1, d, J])/(m^4-u^2))
X[k_,u_,m_,J_,d_]:=(2m^2+u)/(u m^2 (m^2+u)) (m^2 Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2u)/m^2))/2])/(u m^2 (m^2+u))^(k/2)-Residue[((2m^2+t)(m^2-t)(m^2+2t))/(m^2 (u-t)t(m^2-u)(m^2+t)(m^2+u+t)) (m^2 Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-(1+(2t)/m^2))/2])/(t m^2 (m^2+t))^(k/2),{t,0}]

K1[n_, J_, q_] := \[Integral]p^n A[m, J, -p^2, 7]\[DifferentialD]p /.{p -> 1 , m -> q};
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
        
        
       {
       
       Poly[0, 1], Poly[0, 1.1], Poly[0, 1.2], Poly[0, 1.3], Poly[0, 1.4], Poly[0, 1.5], Poly[0, 1.6], Poly[0, 1.7], Poly[0, 1.8], Poly[0, 1.9], Poly[0, 2], 
       Poly[0, 1.05], Poly[0, 1.15], Poly[0, 1.25], Poly[0, 1.35], Poly[0, 1.45], Poly[0, 1.55], Poly[0, 1.65], Poly[0, 1.75], Poly[0, 1.85], Poly[0, 1.95],
       Poly[2, 1], Poly[2, 1.1], Poly[2, 1.2], Poly[2, 1.3], Poly[2, 1.4], Poly[2, 1.5], Poly[2, 1.6], Poly[2, 1.7], Poly[2, 1.8], Poly[2, 1.9], Poly[2, 2],
       Poly[2, 1.05], Poly[2, 1.15], Poly[2, 1.25], Poly[2, 1.35], Poly[2, 1.45], Poly[2, 1.55], Poly[2, 1.65], Poly[2, 1.75], Poly[2, 1.85], Poly[2, 1.95],
       Poly[4, 1], Poly[4, 1.1], Poly[4, 1.2], Poly[4, 1.3], Poly[4, 1.4], Poly[4, 1.5], Poly[4, 1.6], Poly[4, 1.7], Poly[4, 1.8], Poly[4, 1.9], Poly[4, 2],
       Poly[4, 1.05], Poly[4, 1.15], Poly[4, 1.25], Poly[4, 1.35], Poly[4, 1.45], Poly[4, 1.55], Poly[4, 1.65], Poly[4, 1.75], Poly[4, 1.85], Poly[4, 1.95],
       Poly[6, 1], Poly[6, 1.1], Poly[6, 1.2], Poly[6, 1.3], Poly[6, 1.4], Poly[6, 1.5], Poly[6, 1.6], Poly[6, 1.7], Poly[6, 1.8], Poly[6, 1.9], Poly[6, 2],  
       Poly[8, 1], Poly[8, 1.1], Poly[8, 1.2], Poly[8, 1.3], Poly[8, 1.4], Poly[8, 1.5], Poly[8, 1.6], Poly[8, 1.7], Poly[8, 1.8], Poly[8, 1.9], Poly[8, 2], 
       Poly[10, 1], Poly[10, 1.1], Poly[10, 1.2], Poly[10, 1.3], Poly[10, 1.4], Poly[10, 1.5], Poly[10, 1.6], Poly[10, 1.7], Poly[10, 1.8], Poly[10, 1.9], Poly[10, 2], 
       Poly[12, 1], Poly[12, 1.1], Poly[12, 1.2], Poly[12, 1.3], Poly[12, 1.4], Poly[12, 1.5], Poly[12, 1.6], Poly[12, 1.7], Poly[12, 1.8], Poly[12, 1.9], Poly[12, 2],
       Poly[14, 1], Poly[14, 1.1], Poly[14, 1.2], Poly[14, 1.3], Poly[14, 1.4], Poly[14, 1.5], Poly[14, 1.6], Poly[14, 1.7], Poly[14, 1.8], Poly[14, 1.9], Poly[14, 2],
       Poly[16, 1], Poly[16, 1.1], Poly[16, 1.2], Poly[16, 1.3], Poly[16, 1.4], Poly[16, 1.5], Poly[16, 1.6], Poly[16, 1.7], Poly[16, 1.8], Poly[16, 1.9], Poly[16 , 2],
       Poly[18, 1], Poly[18, 2], Poly[20, 1], Poly[20, 2], Poly[22, 1], Poly[22, 2], Poly[24, 1], Poly[24, 2], Poly[26, 1], Poly[26, 2], Poly[28, 1], Poly[28, 2], Poly[30, 1], Poly[30, 2], 
       Poly[32, 1], Poly[32, 2], Poly[34, 1], Poly[34, 2], Poly[36, 1], Poly[36, 2], Poly[38, 1], Poly[38, 2], Poly[40, 1], Poly[40, 2], Poly[400, 1], Poly[400, 2](*,
       NewPoly[10], NewPoly[15], NewPoly[20], NewPoly[25], NewPoly[30], NewPoly[35], NewPoly[40]*)
        (*
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 0, 1], K[3, 0, 1], K[4, 0, 1], K[5, 0, 1], K[6, 0, 1], K[7, 0, 1], M[0, 0, 1], M[1, 0, 1], M[2, 0, 1], P[0, 0, 1], P[1, 0, 1], P[2, 0, 1]}}}],
 
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 0, 1.3], K[3, 0, 1.3], K[4, 0, 1.3], K[5, 0, 1.3], K[6, 0, 1.3], K[7, 0, 1.3], M[0, 0, 1.3], M[1, 0, 1.3], M[2, 0, 1.3], P[0, 0, 1.3], P[1, 0, 1.3], P[2, 0, 1.3]}}}],
                      
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 0, 1.5], K[3, 0, 1.5], K[4, 0, 1.5], K[5, 0, 1.5], K[6, 0, 1.5], K[7, 0, 1.5], M[0, 0, 1.5], M[1, 0, 1.5], M[2, 0, 1.5], P[0, 0, 1.5], P[1, 0, 1.5], P[2, 0, 1.5]}}}],
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 2, 1], K[3, 2, 1], K[4, 2, 1], K[5, 2, 1], K[6, 2, 1], K[7, 2, 1], M[0, 2, 1], M[1, 2, 1], M[2, 2, 1], P[0, 2, 1], P[1, 2, 1], P[2, 2, 1]}}}],        
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 2, 1.5], K[3, 2, 1.5], K[4, 2, 1.5], K[5, 2, 1.5], K[6, 2, 1.5], K[7, 2, 1.5], M[0, 2, 1.5], M[1, 2, 1.5], M[2, 2, 1.5], P[0, 2, 1.5], P[1, 2, 1.5], P[2, 2, 1.5]}}}],
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 4, 1], K[3, 4, 1], K[4, 4, 1], K[5, 4, 1], K[6, 4, 1], K[7, 4, 1], M[0, 4, 1], M[1, 4, 1], M[2, 4, 1], P[0, 4, 1], P[1, 4, 1], P[2, 4, 1]}}}],
 
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 4, 1.5], K[3, 4, 1.5], K[4, 4, 1.5], K[5, 4, 1.5], K[6, 4, 1.5], K[7, 4, 1.5], M[0, 4, 1.5], M[1, 4, 1.5], M[2, 4, 1.5], P[0, 4, 1.5], P[1, 4, 1.5], P[2, 4, 1.5]}}}],
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 6, 1], K[3, 6, 1], K[4, 6, 1], K[5, 6, 1], K[6, 6, 1], K[7, 6, 1], M[0, 6, 1], M[1, 6, 1], M[2, 6, 1], P[0, 6, 1], P[1, 6, 1], P[2, 6, 1]}}}],
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 6, 1.5], K[3, 6, 1.5], K[4, 6, 1.5], K[5, 6, 1.5], K[6, 6, 1.5], K[7, 6, 1.5], M[0, 6, 1.5], M[1, 6, 1.5], M[2, 6, 1.5], P[0, 6, 1.5], P[1, 6, 1.5], P[2, 6, 1.5]}}}],
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 8, 1.5], K[3, 8, 1.5], K[4, 8, 1.5], K[5, 8, 1.5], K[6, 8, 1.5], K[7, 8, 1.5], M[0, 8, 1.5], M[1, 8, 1.5], M[2, 8, 1.5], P[0, 8, 1.5], P[1, 8, 1.5], P[2, 8, 1.5]}}}],       
        
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 8, 1], K[3, 8, 1], K[4, 8, 1], K[5, 8, 1], K[6, 8, 1], K[7, 8, 1], M[0, 8, 1], M[1, 8, 1], M[2, 8, 1], P[0, 8, 1], P[1, 8, 1], P[2, 8, 1]}}}],       
                             
        PositiveMatrixWithPrefactor[
        DampedRational[1,{},1/E,x],{{{K[2, 10, 1], K[3, 10, 1], K[4, 10, 1], K[5, 10, 1], K[6, 10, 1], K[7, 10, 1], M[0, 10, 1], M[1, 10, 1], M[2, 10, 1], P[0, 10, 1], P[1, 10, 1], P[2, 10, 1]}}}]
        *)
        
        },

        
      
        norm =  {1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        obj  = -1 * ({1, 1/2, 1/3, 1/4, 1/5, 1/6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} + 2 * (-2.8)  * {1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
        (*
        norm =  {1/5, 1/6, 1/7, 1/8, 1/9, 1/10},
        obj  =  {1, 1/2, 1/3, 1/4, 1/5, 1/6}
        *)       
    },

    Print[norm];
    Print[obj];
    WritePmpJson[datfile, SDP[obj, norm, pols]]]

TSDP["out.json"]
