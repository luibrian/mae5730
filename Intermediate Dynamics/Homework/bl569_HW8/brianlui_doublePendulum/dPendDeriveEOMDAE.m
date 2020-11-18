%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the double pendulum. This does it using DAE.

syms d1 d2 L1 m1 I1 m2 I2 g real           %parameters (Note d one, not d L)
syms theta1 theta2 theta1dot theta2dot real              %states
syms theta1doubledot theta2doubledot real                %states
syms x1 x1dot x1doubledot y1 y1dot y1doubledot real      %states
syms x2 x2dot x2doubledot y2 y2dot y2doubledot real      %states
syms e1 e2 e3 er1 etheta1 er2 etheta2 real               %unit vectors
syms rG1relO aG1relO real                                %relevant vectors
syms rG2relO aG2relO rG2relE aG2relE rErelO aErelO real  %relevant vectors
syms MrelO MrelE H1relODot H2relODot H real
syms ROx ROy REx REy real                                %Reaction forces

%10 unknowns, 10 equations (6 from LMB/AMB, 4 from constraint)
unknowns = [x1doubledot, y1doubledot, theta1doubledot, ...
    x2doubledot, y2doubledot, theta2doubledot, ...
    ROx, ROy, REx, REy];

%defining unit vectors
e1 = [1 0 0]'; 
e2 = [0 1 0]'; 
e3 = cross(e1,e2); 
er1 = cos(theta1)*e1 + sin(theta1)*e2;
etheta1 = cross(e3, er1);
er2 = cos(theta2)*e1 + sin(theta2)*e2;
etheta2 = cross(e3, er2);

%defining vectors of interest
%For the first bar
rG1relO = d1*er1;
aG1relO = d1*theta1doubledot *etheta1 - d1*theta1dot^2 * er1;

%For the second bar (forearm)
rG2relE = d2*er2;
rErelO = L1*er1;
rG2relO = rG2relE + rErelO;
rErelG2 = -rG2relE;

rOrelG1 = -rG1relO;
rErelG1 = rErelO - rG1relO;

aG2relE = d2*theta2doubledot*etheta2 - d2*theta2dot^2 * er2;
aErelO  = L1*theta1doubledot*etheta1 - L1*theta1dot^2 * er1;
aG2relO = aG2relE + aErelO;

%defining forces acting on system
%don't define reaction forces because we ignore them when we do AMB
gravityForce1 = m1*g*e1;
gravityForce2 = m2*g*e1;

%Reaction forces
reactO = ROx*e1 + ROy*e2;
reactE = REx*e1 + REy*e2;

%Forces on Links
forcesOnFirstLink = -reactO + reactE +gravityForce1;
forcesOnSecondLink = -reactE + gravityForce2;

%Linear Momentum Balance Equations
%First 2 equations (actually more like 4 equations
%lmb1 = forcesOnFirstLink - m1*aG1relO;
lmb1x = forcesOnFirstLink(1) - m1*x1doubledot;   %LMB1 in the x direction
lmb1y = forcesOnFirstLink(2) - m1*y1doubledot;   %LMB1 in the y direction
%lmb2 = forcesOnSecondLink - m2*aG2relO;
lmb2x = forcesOnSecondLink(1) - m2*x2doubledot;   %LMB2 in the x direction
lmb2y = forcesOnSecondLink(2) - m2*y2doubledot;   %LMB2 in the y direction

%Angular Momentum Balance Equations
HrelODot = cross(rG1relO, m1*aG1relO) + I1*theta1doubledot*e3;
MrelO = cross(rG1relO, gravityForce1) + cross(rG2relO, gravityForce2);

HrelG1Dot = I1*theta1doubledot*e3;
MrelG1 = cross(rOrelG1, -reactO) + cross(rErelG1, reactE)

HrelG2Dot = I2*theta2doubledot*e3;
MrelG2 = cross(rErelG2, -reactE);

%amb1 = dot(HrelODot - MrelO, e3)
amb1 = dot(HrelG1Dot - MrelG1, e3)
amb2 = dot(HrelG2Dot - MrelG2, e3)

%Constraint Equations
constraint1x = x1doubledot - aG1relO(1);
constraint1y = y1doubledot - aG1relO(2);

constraint2x = x2doubledot - aG2relO(1);
constraint2y = y2doubledot - aG2relO(2);



%Putting all the equations into a matrix
eqns = [lmb1x, lmb1y lmb2x, lmb2y amb1, amb2,...
    constraint1x, constraint1y, constraint2x, constraint2y];

[A,b] = equationsToMatrix(eqns, unknowns);
A = simplify(A)
size(A)
b = simplify(b)

matlabFunction(A, 'file', 'rhsStuffmassMatrixDAE');
matlabFunction(b, 'file', 'rhsStuffbVectorDAE');