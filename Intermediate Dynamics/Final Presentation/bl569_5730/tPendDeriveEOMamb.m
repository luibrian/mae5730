%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the double pendulum. This does it using AMB. 

syms d1 d2 d3 L1 L2 m1 I1 m2 I2 m3 I3 g real     %parameters (Note d one, not d L)
syms theta1 theta2 theta3 theta1dot theta2dot theta3dot real   %states
syms theta1doubledot theta2doubledot theta3doubledot real      %states
syms e1 e2 e3 er1 etheta1 er2 etheta2 er3 etheta3 real   %unit vectors
syms rG1relO aG1relO                                %relevant vectors
syms rG2relO aG2relO rG2relE aG2relE rErelO aErelO  %relevant vectors
syms rG3relO aG3relO rG3relF aG3relF rFrelE aFrelE real  %relveant vectors 
syms rG3relE aG3relE real
syms MrelO MrelE MrelF H1relODot H2relODot H3relODot real

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
aG1relO = d1*theta1doubledot *etheta1 - d1*theta1dot^2 * er1;

%For the second bar (forearm)
rG2relE = d2*er2;
rErelO = L1*er1;
rG2relO = rG2relE + rErelO;

aG2relE = d2*theta2doubledot*etheta2 - d2*theta2dot^2 * er2;
aErelO  = L1*theta1doubledot*etheta1 - L1*theta1dot^2 * er1;
aG2relO = aG2relE + aErelO;

%For the third bar
rG3relF = d3*er3; 
rFrelG3 = -rG3relF;

rFrelE  = L2*er2;
rFrelO = rFrelE + rErelO;
    
rG3relE = rG3relF + rFrelE;
rG3relO = rG3relF + rFrelE + rErelO;

rFrelG2 = rFrelO - rG2relO;

aG3relF = d3*theta3doubledot*etheta3 - d3*theta3dot^2 * er3;
aFrelE = L2*theta2doubledot*etheta2 - L2*theta2dot^2 * er2;
aG3relE = aG3relF + aFrelE;
aG3relO = aG3relF + aFrelE + aErelO;



%defining forces acting on system
%don't define reaction forces because we ignore them when we do AMB
gravityForce1 = m1*g*e1;
gravityForce2 = m2*g*e1;
gravityForce3 = m3*g*e1;

%The stuff that goes into the angular momentum balance
%AMB about point O
H1relODot = cross(rG1relO, m1*aG1relO) + I1*theta1doubledot*e3;
H2relODot = cross(rG2relO, m2*aG2relO) + I2*theta2doubledot*e3;
H3relODot = cross(rG3relO, m3*aG3relO) + I3*theta3doubledot*e3;
HrelODot = H1relODot + H2relODot + H3relODot;

%AMB about point E
H2relEDot = cross(rG2relE, m2*aG2relO) + I2*theta2doubledot*e3;
H3relEDot = cross(rG3relE, m3*aG3relO) + I3*theta3doubledot*e3;
HrelEDot = H2relEDot + H3relEDot;

%AMB about point F
HrelFDot = cross(rG3relF, m2*aG3relO) + I3*theta3doubledot*e3;

%The net moments of each of these cases
MrelO = cross(rG1relO, gravityForce1) + cross(rG2relO, gravityForce2) +...
    cross(rG3relO, gravityForce3);
MrelE = cross(rG2relE, gravityForce2) + cross(rG3relE, gravityForce3);
MrelF = cross(rG3relF, gravityForce3);


%equations to solve
%These equations should equal 0 at all times for system
%eqn1 = dot(MrelO - HrelODot, e3)    %I think this also works?
eqn1 = MrelO - HrelODot;
eqn1 = eqn1(3);
eqn2 = MrelE - HrelEDot;
eqn2 = eqn2(3);
eqn3 = MrelF - HrelFDot;
eqn3 = eqn3(3);
eqns = [eqn1; eqn2; eqn3];

thetaDoubleDots = [theta1doubledot, theta2doubledot, theta3doubledot];
[A,b] = equationsToMatrix(eqns, thetaDoubleDots);
A = simplify(A)
b = simplify(b)


matlabFunction(A, 'file', 'rhsStuffmassMatrixAMB');
matlabFunction(b, 'file', 'rhsStuffbVectorAMB');