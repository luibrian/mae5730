%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the triple pendulum. This does it using DAE.

syms d1 d2 d3 L1 L2 L3 m1 I1 m2 I2 m3 I3 g real     %parameters (Note d one, not d L)
syms theta1 theta2 theta3 theta1dot theta2dot theta3dot real   %states
syms theta1doubledot theta2doubledot theta3doubledot real      %states
syms x1 x1dot x1doubledot y1 y1dot y1doubledot real      %states
syms x2 x2dot x2doubledot y2 y2dot y2doubledot real      %states
syms x3 x3dot x3doubledot y3 y3dot y3doubledot real      %states
syms e1 e2 e3 er1 etheta1 er2 etheta2 er3 etheta3 real   %unit vectors
syms rG1relO aG1relO real                                %relevant vectors
syms rG2relO aG2relO rG2relE aG2relE rErelO aErelO real  %relevant vectors
syms rG3relO aG3relO rG3relF aG3relF rFrelE aFrelE real  %relveant vectors 
syms MrelG1 MrelG2 MrelG3 real
syms ROx ROy REx REy RFx RFy real                        %Reaction forces

%15 unknowns, 15 equations (9 from LMB/AMB, 6 from constraint)
unknowns = [x1doubledot, y1doubledot, theta1doubledot, ...
    x2doubledot, y2doubledot, theta2doubledot, ...
    x3doubledot, y3doubledot, theta3doubledot, ...
    ROx, ROy, REx, REy, RFx, RFy];

%defining unit vectors
e1 = [1 0 0]'; 
e2 = [0 1 0]'; 
e3 = cross(e1,e2); 

er1 = cos(theta1)*e1 + sin(theta1)*e2;
etheta1 = cross(e3, er1);
er2 = cos(theta2)*e1 + sin(theta2)*e2;
etheta2 = cross(e3, er2);
er3 = cos(theta3)*e1 + sin(theta3)*e2;
etheta3 = cross(e3, er3);

%defining vectors of interest
%For the first bar
rG1relO = d1*er1;
rOrelG1 = -rG1relO;
aG1relO = d1*theta1doubledot *etheta1 - d1*theta1dot^2 * er1;

%For the second bar (forearm)
rG2relE = d2*er2;
rErelG2 = -rG2relE;
rErelO = L1*er1;
rG2relO = rG2relE + rErelO;

aG2relE = d2*theta2doubledot*etheta2 - d2*theta2dot^2 * er2;
aErelO  = L1*theta1doubledot*etheta1 - L1*theta1dot^2 * er1;
aG2relO = aG2relE + aErelO;

rErelG1 = rErelO - rG1relO;
%aErelG1 = aErelO - aG1relO;     %Added this

%For the third bar
rG3relF = d3*er3; 
rFrelG3 = -rG3relF;

rFrelE  = L2*er2;
rFrelO = rFrelE + rErelO;

rG3relE = rG3relF + rFrelE;
rG3relO = rG3relF + rFrelE + rErelO;

rFrelG2 = rFrelO - rG2relO;
%added these here too
%rHrelF = L3*er3;
%rHrelO = rHrelF + rFrelO;
%rHrelG3 = rHrelO - rG3relO;

aG3relF = d3*theta3doubledot*etheta3 - d3*theta3dot^2 * er3;
%aG3relF = cross(theta3doubledot*e3, rG3relF) - theta3dot^2*d3*er3
aFrelE = L2*theta2doubledot*etheta2 - L2*theta2dot^2 * er2;
aG3relE = aG3relF + aFrelE;
aG3relO = aG3relF + aFrelE + aErelO;
%aFrelG2 = aFrelE - aG2relE;         %added this
%aHrelF = L3*theta3doubledot*etheta3 - L3*theta3dot^2 * er3;     %added this
%aHrelG3 = aHrelF - aG3relF;         %added this
 

%defining forces acting on system
gravityForce1 = m1*g*e1;
gravityForce2 = m2*g*e1;
gravityForce3 = m3*g*e1;

%Reaction forces
reactO = ROx*e1 + ROy*e2;
reactE = REx*e1 + REy*e2;
reactF = RFx*e1 + RFy*e2;

%Forces on Links
forcesOnFirstLink  = -reactO + reactE + gravityForce1;
forcesOnSecondLink = -reactE + reactF + gravityForce2;
forcesOnThirdLink  = -reactF + gravityForce3;

%Linear Momentum Balance Equations
%First 2 equations (actually more like 6 equations)
%lmb1 = forcesOnFirstLink - m1*aG1relO;
lmb1x = forcesOnFirstLink(1) - m1*x1doubledot;   %LMB1 in the x direction
lmb1y = forcesOnFirstLink(2) - m1*y1doubledot;   %LMB1 in the y direction
%lmb2 = forcesOnSecondLink - m2*aG2relO;
lmb2x = forcesOnSecondLink(1) - m2*x2doubledot;   %LMB2 in the x direction
lmb2y = forcesOnSecondLink(2) - m2*y2doubledot;   %LMB2 in the y direction
%lmb3 = forcesonThirdLink - m3*aG3relO;
lmb3x = forcesOnThirdLink(1) - m3*x3doubledot;
lmb3y = forcesOnThirdLink(2) - m3*y3doubledot;

%Angular Momentum Balance Equations
%HrelODot = cross(rG1relO, m1*aG1relO) + I1*theta1doubledot*e3;
%MrelO = cross(rG1relO, gravityForce1) + cross(rG2relO, gravityForce2);
%
HrelG1Dot = I1*theta1doubledot*e3;
MrelG1 = cross(rOrelG1, -reactO) + cross(rErelG1, reactE);

HrelG2Dot = I2*theta2doubledot*e3;
MrelG2 = cross(rErelG2, -reactE) + cross(rFrelG2, reactF);

HrelG3Dot = I3*theta3doubledot*e3;
MrelG3 = cross(rFrelG3, -reactF);

%amb1 = dot(HrelODot - MrelO, e3)
amb1 = dot(HrelG1Dot - MrelG1, e3);
amb2 = dot(HrelG2Dot - MrelG2, e3);
amb3 = dot(HrelG3Dot - MrelG3, e3);
%

%{
%The stuff that goes into the angular momentum balance
%AMB about point O, ie. these equations also work.
HrelODot = cross(rG1relO, m1*aG1relO) + I1*theta1doubledot*e3;
%AMB about point E
HrelEDot = cross(rG2relE, m2*aG2relO) + I2*theta2doubledot*e3;
%AMB about point F
HrelFDot = cross(rG3relF, m2*aG3relO) + I3*theta3doubledot*e3;

MrelO = cross(rG1relO, gravityForce1) + cross(rErelO, reactE);
MrelE = cross(rG2relE, gravityForce2) + cross(rFrelE, reactF);
MrelF = cross(rG3relF, gravityForce3);

amb1 = dot(HrelODot - MrelO, e3);
amb2 = dot(HrelEDot - MrelE, e3);
amb3 = dot(HrelFDot - MrelF, e3);
%}

%Constraint Equations
constraint1x = x1doubledot - aG1relO(1);
constraint1y = y1doubledot - aG1relO(2);

constraint2x = x2doubledot - aG2relO(1);
constraint2y = y2doubledot - aG2relO(2);

constraint3x = x3doubledot - aG3relO(1);
constraint3y = y3doubledot - aG3relO(2);

%{
%These also work
constraint1x = x1doubledot - aG1relO(1)
constraint1y = y1doubledot - aG1relO(2);

constraint2x = x2doubledot - aG2relE(1) - (x1doubledot + aErelG1(1));
constraint2y = y2doubledot - aG2relE(2) - (y1doubledot + aErelG1(2));

constraint3x = x3doubledot - aG3relF(1) - (x2doubledot + aFrelG2(1))
constraint3y = y3doubledot - aG3relF(2) - (y2doubledot + aFrelG2(2));
%}

%Putting all the equations into a matrix 15 equations
eqns = [lmb1x, lmb1y, lmb2x, lmb2y, lmb3x, lmb3y,...
    amb1, amb2, amb3,...
    constraint1x, constraint1y,...
    constraint2x, constraint2y,...
    constraint3x, constraint3y];

[A,b] = equationsToMatrix(eqns, unknowns);
A = simplify(A)
b = simplify(b)

matlabFunction(A, 'file', 'rhsStuffmassMatrixDAE');
matlabFunction(b, 'file', 'rhsStuffbVectorDAE');