% six bar linkage
% Static Equilibrium

clc;
clear;

% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

% define the lengths of the bars

lAB = norm(B - A);
lBC = norm(C - B);
lCE = norm(C - E);
lCD = norm(D - C);
lBE = norm(E - B);
lEF = norm(F - E);
lFG = norm(G -F);

% Weight od each link
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

%Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];

InputTorque = [0 0 Tin];

% Applied Force
AppliedForce = [50 0 0];

% Static Equilibrium conditions for Link AB
%Fa + Fb + WeightofAB
eqn1 = ForceA +ForceB + WAB == 0;

%Sum of Moments = 0 with respect to CoM of link AB
% S1A x FA + S1B x FB + InputTorque = 0

eqn2 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == 0;

%equations for Link BEC
%Sum of forces =0
%-Fb +Fc +Fe +WBEC = 0

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

%Sum of Moments =0
% S2B x FB + S2C x FC +S2E x FE = 0
eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;

% Link CD
eqn5 = -ForceC +ForceD +WCD == 0;

% Sum of moments
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;

% Link EF
eqn7 = -ForceE +ForceF + WEF == 0;

% Moments
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

% Link FG
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

% Moments
eqn10 = cross(F-S5, -ForceF) + cross(G-S2, ForceG) == 0;

% Solving the 10 equations

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);

% Angular Velocity Calculations
% Loop ABCDA


syms wBEC wCD 
omega_AB =  [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB, B-A) + cross(omega_BEC, C-B) +cross(omega_CD, D-C) == 0;
loop1Solution = solve(eqn11,[wBEC wCD]);

% Extraxt angular velocities from the loop solution
angularVelocity_BEC = double(loop1Solution.wBEC);
angularVelocity_CD = double(loop1Solution.wCD);

omegaBEC = [0 0 angularVelocity_BEC];
omegaCD = [0 0 angularVelocity_CD];

% Loop DCEFGD

syms wEF wFG
omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG]

eqn12 = cross(omegaCD, C-D) + cross(omegaBEC, E-C) + cross(omega_EF, F-E) + cross(omega_FG, G-F) == 0;

loop2Solution = solve(eqn12, [wEF wFG]);

angularVelocity_EF = double(loop2Solution.wEF);
angularVelocity_FG = double(loop2Solution.wFG);

% Extraxt angular accelerations from the loop solution
% Loop ABCDA

syms aBEC aCD aAB
alpha_AB = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];

a_C_B = cross(alpha_BEC, C-B) + cross(omegaBEC, cross(omegaBEC, C-B));
a_D_C = cross(alpha_CD, D-C) + cross(omegaCD, cross(omegaCD, D-C));

eqn13 = a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13, [aBEC aCD]);

alphaBEC = double(loop1AccSolution.aBEC);
alphaCD = double(loop1AccSolution.aCD);

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

% Loop DCEFGD
syms aEF aFG
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

% a_C_D + a_E_C + a_F_E + a_G_F = 0
a_C_D = cross(alphaCD_vector, C-D) + cross(omegaCD, cross(omegaCD, C-D));
a_E_C = cross(alphaBEC_vector, E-C) + cross(omegaBEC, cross(omegaBEC, E-C));

angVel_EF = [0 0 angularVelocity_EF];
angVel_FG = [0 0 angularVelocity_FG];

a_F_E = cross(alpha_EF, F-E) + cross(angVel_EF, cross(angVel_EF, F-E));
a_G_F = cross(alpha_FG, G-F) + cross(angVel_FG, cross(angVel_FG, G-F));

eqn13 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn13, [aEF aFG]);

alphaEF = double(loop2AccSolution.aEF)
alphaFG = double(loop2AccSolution.aFG)

% Velocity at Joint S4

vB_A = cross(omega_AB, B-A);
vE_B = cross(omegaBEC, E-B);

% vE_A = VEB + VBA

vE_A = vE_B +vB_A;

%V_S4/G = V_S4_F +V_F_G
V_S4_F = cross(angVel_EF, S4-F);
V_F_G = cross(angVel_FG, F-G);

vS4_G = V_S4_F +V_F_G;

% Velocity at S1
V_S1_A = cross(WAB, S1-A);

% Velocity at S2

% VS2_A = VS2_B + VB_A
V_S2_B = cross(omegaBEC, S2-B);

V_S2_A = V_S2_B + vB_A;

% Velocity at S3
V_S3_D = cross(omegaCD, S3-D);

% Velocity of S5
V_S5_G = cross(angVel_FG, S5-G);

% Calculating Linear Accelerations

% Acceleration at S1
A_S1_A = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));

% Acceleration at S2
% A_S2A = A_S2B + A_BA and A_BA = 0 since angular velocity is constant

A_S2_A = cross(alphaBEC_vector, S2-B) + cross(omegaBEC, cross(omegaBEC, S2-B));

% Acceleration at S3
A_S3_D = cross(alphaCD_vector, S3-D) + cross(omegaCD, cross(omegaCD, S3-D));

% Acceleration at S4
% A_S4G = A_S4F + alphaFG
A_S4_F = cross([0 0 alphaEF], F-S4) + cross(omega_EF, cross(omega_EF, F-S4));
A_S4_G = A_S4_F + alphaFG;

% Acceleration at S5
A_S5_G = cross([0 0 alphaFG], G-S5) + cross(omega_FG, cross(omega_FG, G-S5));

%Newton's Second Moment 
MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;

% Mass Moments of Inertia
J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

% Define Forces

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NInputTorque = [0 0 NTin];

%Equations for link AB
% Sum of Forces
eqn15 = NForceA + NForceB + WAB == MassAB * A_S1_A;

% Sum of Moments
eqn16 = cross(A-S1, NForceA) + cross(B-S1, NForceB) + NInputTorque == J_AB * A_S1_A;

% Link BEC
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC * A_S2_A;
eqn18 = cross(B-S2, -NForceB) + cross(C-S2, NForceC) + cross(E-S2, NForceE) == J_BEC * alphaBEC_vector;

% Link CD
eqn19 = -NForceC + NForceD +WCD == MassCD * A_S3_D;
eqn20 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) == J_CD * alphaCD_vector;

% Link EF
eqn21 = -NForceE + NForceF + WEF == MassEF * A_S4_G;
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF * alphaEF * [0 0 alphaEF];

eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * A_S5_G;
eqn24 = cross(F-S5, -NForceF) + cross(G-S5, NForceG) == J_FG * [0 0 alphaFG];

NeqnMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];
DynamicSolution = solve(NeqnMatrix, [NFAx, NFAy, NFBx, NFBy, NFCx, NFCy, NFDx, NFDy, NFEx, NFEy, NFFx, NFFy, NFGx, NFGy, NTin]);

NForce_Ax = double(DynamicSolution.NFAx);
NForce_Ay = double(DynamicSolution.NFAy);
NForce_Bx = double(DynamicSolution.NFBx);
NForce_By = double(DynamicSolution.NFBy);
NForce_Cx = double(DynamicSolution.NFCx);
NForce_Cy = double(DynamicSolution.NFCy);
NForce_Dx = double(DynamicSolution.NFDx);
NForce_Dy = double(DynamicSolution.NFDy);
NForce_Ex = double(DynamicSolution.NFEx);
NForce_Ey = double(DynamicSolution.NFEy);
NForce_Fx = double(DynamicSolution.NFFx);
NForce_Fy = double(DynamicSolution.NFFy);
NForce_Gx = double(DynamicSolution.NFGx);
NForce_Gy = double(DynamicSolution.NFGy);