%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the double pendulum. This does it using AMB. 

syms d1 d2 L1 m1 I1 m2 I2 g           %parameters (Note d one, not d L)
syms theta1 theta2 theta1dot theta2dot theta1doubledot theta2doubledot %states
syms e1 e2 e3 er1 etheta1 er2 etheta2               %unit vectors
syms rG1relO aG1relO                                %relevant vectors
syms rG2relO aG2relO rG2relE aG2relE rErelO aErelO  %relevant vectors
syms MrelO MrelE H1relODot H2relODot H

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

aG2relE = d2*theta2doubledot*etheta2 - d2*theta2dot^2 * er2;
aErelO  = L1*theta1doubledot*etheta1 - L1*theta1dot^2 * er1;
aG2relO = aG2relE + aErelO;

%defining forces acting on system
%don't define reaction forces because we ignore them when we do AMB
gravityForce1 = m1*g*e1;
gravityForce2 = m2*g*e1;

%The stuff that goes into the angular momentum balance
H1relODot = cross(rG1relO, m1*aG1relO) + I1*theta1doubledot*e3;
H2relODot = cross(rG2relO, m2*aG2relO) + I2*theta2doubledot*e3;
HrelODot = H1relODot + H2relODot;

HrelEDot = cross(rG2relE, m2*aG2relO) + I2*theta2doubledot*e3;

MrelO = cross(rG1relO, gravityForce1) + cross(rG2relO, gravityForce2);
MrelE = cross(rG2relE, gravityForce2);


%equations to solve
%These equations should equal 0 at all times for system
%eqn1 = dot(MrelO - HrelODot, e3)    %I think this also works?
eqn1 = MrelO - HrelODot;
eqn1 = eqn1(3);
eqn2 = MrelE - HrelEDot;
eqn2 = eqn2(3);
eqns = [eqn1; eqn2];

% %Option 1: using solve
% [alpha1, alpha2] = solve([eqn1, eqn2], [theta1doubledot, theta2doubledot]);
% alpha1 = simplify(alpha1);
% alpha2 = simplify(alpha2);
% alphas = [alpha1,;alpha2];
% 
% matlabFunction(alphas, 'file', 'rhsStuffAMB');


%Option 2: Using jacobian
%thetaDoubleDots = [theta1doubledot, theta2doubledot];
%jac = jacobian(eqns, thetaDoubleDots);
%A = simplify(jac);
%b = simplify(A*thetaDoubleDots' - eqns);

%Option 3: Use equationsToMatrix
thetaDoubleDots = [theta1doubledot, theta2doubledot];
[A,b] = equationsToMatrix(eqns, thetaDoubleDots);
A = simplify(A)
b = simplify(b)


matlabFunction(A, 'file', 'rhsStuffmassMatrixAMB');
matlabFunction(b, 'file', 'rhsStuffbVectorAMB');