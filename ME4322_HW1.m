clc;
clear;

% define joints
A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];

%Where force is applied
P = [-1.679 4.403 0];

% joint values
new_B_x(1) = B(1);
new_B_y(1) = B(1);
new_C_x(1) = C(1);
new_C_y(1) = C(1);
new_E_x(1) = E(1);
new_E_y(1) = E(1);
new_F_x(1) = F(1);
new_F_y(1) = F(1);

% lengths of links
lAB = norm(B-A);
lBC = norm(C-B);
lCD = norm(D-C);
lDE = norm(E-D);
lEF = norm(F-E);
lFG = norm(G-F);

% Weights of links
WAB= [0 -228.942837 0];
WBC= [0 -553.5722178 0];
WDE= [0 -965.914182 0];
WEF= [0 -461.046456 0];
WFG= [0 -1003.455875 0];

%Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C)/2;
S3 = (D+E)/2;
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
AppliedForce = [0 -200 0];

% Static Equilibrium conditions for Link AB
%Fa + Fb + WeightofAB
eqn1 = ForceA +ForceB + WAB == 0;

%Sum of Moments = 0 with respect to CoM of link AB
% S1A x FA + S1B x FB + InputTorque = 0

eqn2 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == 0;

%equations for Link BC
%Sum of forces =0
%-Fb +Fc +WBEC = 0

eqn3 = -ForceB + ForceC + WBC == 0;

%Sum of Moments =0
% S2B x FB + S2C x FC +S2E x FE = 0
eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) == 0;

% Link DE
eqn5 = ForceD + ForceE - ForceC + WDE == 0;

% Sum of moments
eqn6 = cross(D-S3, ForceD) + cross(E-S3, ForceE) + cross(C-S3, ForceC) == 0;

% Link EF
eqn7 = -ForceE +ForceF + WEF == 0;

% Moments
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

% Link FG
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

% Moments
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) + cross(P-S5, AppliedForce) == 0;

% Solving the 10 equations

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax(1) = double(StaticSolution.FAx);
Force_Ay(1) = double(StaticSolution.FAy);
Force_Bx(1) = double(StaticSolution.FBx);
Force_By(1) = double(StaticSolution.FBy);
Force_Cx(1) = double(StaticSolution.FCx);
Force_Cy(1) = double(StaticSolution.FCy);
Force_Dx(1) = double(StaticSolution.FDx);
Force_Dy(1) = double(StaticSolution.FDy);
Force_Ex(1) = double(StaticSolution.FEx);
Force_Ey(1) = double(StaticSolution.FEy);
Force_Fx(1) = double(StaticSolution.FFx);
Force_Fy(1) = double(StaticSolution.FFy);
Force_Gx(1) = double(StaticSolution.FGx);
Force_Gy(1) = double(StaticSolution.FGy);
TInput(1) = double(StaticSolution.Tin);

% Angular Velocities

syms wBC wDE wEF wFG

% HELP WHERE TO FIND omegaAB
omega_AB = [0 0 2.424];
omega_BC = [0 0 wBC];
omega_DE = [0 0 wDE];
omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];

% Loop ABCD

eqn11 = cross(omega_AB, B-A) + cross(omega_BC, C-B) +cross(omega_DE, E-D) == 0;
loop1Solution = solve(eqn11,[wBC wDE]);

% Extract angular velocities from the loop solution
angularVelocity_BC = double(loop1Solution.wBC);
angularVelocity_DE = double(loop1Solution.wDE);

omegaBC = [0 0 angularVelocity_BC];
omegaDE = [0 0 angularVelocity_DE];

% Loop ABCEFGA

eqn12 = cross(omega_AB, B-A) + cross(omegaBC, C-B) + cross(omegaDE, E-C) + cross(omega_EF, F-E) + cross(omega_FG, G-F) == 0;
loop2Solution = solve(eqn12, [wEF wFG]);

angularVelocity_EF = double(loop2Solution.wEF);
angularVelocity_FG = double(loop2Solution.wFG);

omegaEF = [0 0 angularVelocity_EF];
omegaFG = [0 0 angularVelocity_FG];

% Extract angular accelerations from the loop solution
% Loop ABCDA

syms aBC aDE aEF aFG

alpha_AB = [0 0 0];
alpha_BC = [0 0 aBC];
alpha_DE = [0 0 aDE];
alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

a_B_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
a_C_B = cross(alpha_BC, C-B) + cross(omegaBC, cross(omegaBC, C-B));
a_E_D = cross(alpha_DE, E-D) + cross(omegaDE, cross(omegaDE, E-D));

eqn13 = a_B_A + a_C_B + a_E_D == 0;

loop1AccSolution = solve(eqn13, [aBC aDE]);

alphaBC = double(loop1AccSolution.aBC);
alphaDE = double(loop1AccSolution.aDE);

alphaBC_vector = [0 0 alphaBC];
alphaDE_vector = [0 0 alphaDE];

% Loop ABCEFGA

% a_E_D + a_F_E + a_G_F = 0
a_E_D = cross(alphaDE_vector, E-D) + cross(omegaDE, cross(omegaDE, E-D));
a_F_E = cross(alpha_EF, F-E) + cross(omegaEF, cross(omegaEF, F-E));
a_G_F = cross(alpha_FG, G-F) + cross(omegaFG, cross(omegaFG, G-F));

%eqn14 = alphaBC_vector + a_E_C + a_F_E + a_G_F == 0;
eqn14 = a_E_D + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14, [aEF aFG]);

alphaEF = double(loop2AccSolution.aEF);
alphaFG = double(loop2AccSolution.aFG);

alphaEF_vector = [0 0 alphaEF];
alphaFG_vector = [0 0 alphaFG];

% Velocity at Joint S1
vS1_A = cross(omega_AB, S1-A);

% Velocity at S2
% VS2_A = VS2_B + VB_A
vS2_B = cross(omegaBC, S2-B);
vS2_A = vS2_B + cross(omega_AB, B-A);

% Velocity at S3
vS3_D = cross(omegaDE, S3-D);

% Velocity at S4
% VS4_D = VS4_E + VE_D
vS4_E = cross(omegaEF, S4-E);
vS4_D = vS4_E + cross(omegaDE, E-D);

% Velocity at S5
vS5_G = cross(omegaFG, S5-G);


% Calculating Linear Accelerations

% Acceleration at S1 (AB)
aS1_A = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));

% Acceleration at S2 (BC)
a_BA = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
aS2_A = cross(alphaBC_vector, S2-B) + cross(omegaBC, cross(omegaBC, S2-B)) + a_BA;

% Acceleration at S3 (DE)
aS3_D = cross(alphaDE_vector, S3-D) + cross(omegaDE, cross(omegaDE, S3-D));

% Acceleration at S4 (EF)
% A_S4D = A_S4E + alphaDE
aS4_E = cross(alphaEF_vector, S4-E) + cross(omegaEF, cross(omegaEF, S4-E));
aS4_D = aS4_E +  cross(alphaDE_vector, E-D) + cross(omegaDE, cross(omegaDE, E-D));

% Acceleration at S5 (FG)
aS5_G = cross(alphaFG_vector, S5-G) + cross(omegaFG, cross(omegaFG, S5-G));

%Newton's Second Moment 
% Weights of links
MassAB = 23.33770;
MassBC = 56.42938;
MassDE = 98.4622;
MassEF = 46.9976;
MassFG = 102.28908;

% Mass Moments of Inertia
J_AB = 7.28559828;
J_BC = 98.18793512;
J_DE = 518.082331;
J_EF = 56.97713652;
J_FG = 580.7263835;

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

eqn15 = NForceA + NForceB + WAB == MassAB .* aS1_A;

% Sum of Moments
eqn16 = cross(A-S1, NForceA) + cross(B-S1, NForceB) + NInputTorque == J_AB .* aS1_A;

% Link BC
eqn17 = -NForceB + NForceC + WBC == MassBC .* aS2_A;
eqn18 = cross(B-S2, -NForceB) + cross(C-S2, NForceC) + cross(E-S2, NForceE) == J_BC .* alphaBC_vector;

% Link CD
eqn19 = NForceD + NForceE + WDE == MassDE .* aS3_D;
eqn20 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) == J_DE .* alphaDE_vector;

% Link EF
eqn21 = -NForceE + NForceF + WEF == MassEF .* aS4_D;
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF .* alphaEF;

eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG .* aS5_G;
eqn24 = cross(F-S5, -NForceF) + cross(G-S5, NForceG) + cross(P-S5, AppliedForce) == J_FG .* [0 0 alphaFG];


% eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];
% StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

% Force_Ax = double(StaticSolution.FAx);

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

% Circle intersection technique

% Joint coordinates have been defined
% length of links have been defined

% compute initial angle of input link AB

initial_theta = atan2(B(2) - A(2), B(1) - A(1));

if (initial_theta < 0)
    inputAngle = 2*pi +initial_theta;
else
    inputAngle = initial_theta;
end

for theta = 1 : 1 : 360
    % new position of joint B

    B_new = A + [lAB*cos(inputAngle+deg2rad(theta)) lAB*sin(inputAngle+deg2rad(theta)) 0];

    % new position of C
    % with B_new as center, BC as radius
    % with D as center and DC as radius

    [Cx, Cy] = circcirc(B_new(1), B_new(2), lBC, D(1), D(2), lCD);

    % checking if there is a NaN

    circIntersect_x = any(isnan(vpa(Cx)));
    circIntersect_y = any(isnan(vpa(Cy)));

    if circIntersect_x==0 && circIntersect_y==0
        C_1 = [Cx(1) Cy(1) 0];
        C_2 = [Cx(2) Cy(2) 0];

        %distance to determine if C1 or C2 is correct
        dist1 = norm(C_1 - C);
        dist2 = norm(C_2 - C);

        if (dist1<dist2)
        C_new = vpa(C_1);
        else
        C_new = vpa(C_2);
        end

        %solving for Ex, Ey with similar triangles
        % Ex_new = (Cx-Dx)*lDE/lCD
        Ex_new = abs(D(1) - C(1)) * lDE/lCD;
        Ey_new = abs(D(2) - C(2)) * lDE/lCD;
        E_new = [Ex_new Ey_new 0];

        % Circle intersection for F
        [Fx, Fy] = circcirc(Ex_new, Ey_new, lEF, G(1), G(2), lFG);

        circIntersect_x_F = any(isnan(vpa(Fx)));
        circIntersect_y_F = any(isnan(vpa(Fy)));

        if circIntersect_x_F==0 && circIntersect_y_F==0
            F_1 = [Fx(1) Fy(1) 0];
            F_2 = [Fx(2) Fy(2) 0];

            dist1 = norm(F_1 - F);
            dist2 = norm(F_2 - F);

            if (dist1< dist2)
                F_new = vpa(F_1);
            else 
                F_new = vpa(F_2);
            end

            % Store values for plotting

            new_B_x(theta +1) = B_new(1);
            new_B_y(theta +1) = B_new(1);
            new_C_x(theta +1) = C_new(1);
            new_C_y(theta +1) = C_new(1);
            new_E_x(theta +1) = E_new(1);
            new_E_y(theta +1) = E_new(1);
            new_F_x(theta +1) = F_new(1);
            new_F_y(theta +1) = F_new(1);

        else
            fprintf('New position cannot be determined at angle: %d degree', theta);
        end
    else
        fprintf('New position cannot be determined at angle: %d degree', theta);
    end

    eqn1 = ForceA +ForceB + WAB == 0;

    %Sum of Moments = 0 with respect to CoM of link AB
    % S1A x FA + S1B x FB + InputTorque = 0

    eqn2 = cross(A-S1, ForceA) + cross(B_new-S1, ForceB) + InputTorque == 0;

    %equations for Link BEC
    %Sum of forces =0
    %-Fb +Fc +Fe +WBEC = 0

    eqn3 = -ForceB + ForceC + WBC == 0;

    %Sum of Moments =0
    % S2B x FB + S2C x FC +S2E x FE = 0
    eqn4 = cross(B_new-S2, -ForceB) + cross(C_new-S2, ForceC) == 0;

    % Link CD
    eqn5 = -ForceC +ForceD +WDE == 0;

    % Sum of moments
    eqn6 = cross(C_new-S3, -ForceC) + cross(D-S3, ForceD) == 0;

    % Link EF
    eqn7 = -ForceE +ForceF + WEF == 0;

    % Moments
    eqn8 = cross(E_new-S4, -ForceE) + cross(F_new-S4, ForceF) == 0;

    % Link FG
    eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

    % Moments
    eqn10 = cross(F_new-S5, -ForceF) + cross(G-S2, ForceG) + cross(P-S5, AppliedForce) == 0;

    % Solving the 10 equations

    eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

    StaticSolution = solve(eqnMatrix, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

    Force_Ax(theta + 1) = double(StaticSolution.FAx);
    Force_Ay(theta + 1) = double(StaticSolution.FAy);
    Force_Bx(theta + 1) = double(StaticSolution.FBx);
    Force_By(theta + 1) = double(StaticSolution.FBy);
    Force_Cx(theta + 1) = double(StaticSolution.FCx);
    Force_Cy(theta + 1) = double(StaticSolution.FCy);
    Force_Dx(theta + 1) = double(StaticSolution.FDx);
    Force_Dy(theta + 1) = double(StaticSolution.FDy);
    Force_Ex(theta + 1) = double(StaticSolution.FEx);
    Force_Ey(theta + 1) = double(StaticSolution.FEy);
    Force_Fx(theta + 1) = double(StaticSolution.FFx);
    Force_Fy(theta + 1) = double(StaticSolution.FFy);
    Force_Gx(theta + 1) = double(StaticSolution.FGx);
    Force_Gy(theta + 1) = double(StaticSolution.FGy);
 
end