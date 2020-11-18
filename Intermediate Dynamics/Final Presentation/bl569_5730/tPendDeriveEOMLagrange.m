%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the double pendulum. This does it using Lagrange 
%equations.

%Need the real or else Matlab does some funky stuff with complex conjugates
syms d1 d2 d3 L1 L2 m1 I1 m2 I2 m3 I3 g t real     %parameters (Note d one, not d L)
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

%defining position and velocity vectors of interest
rG1relO = d1*er1;                   %first link

%second link
rErelO  = L1*er1;  
rG2relE = d2*er2;
rG2relO = rErelO + rG2relE;

%third link
rFrelE = L2*er2;
rG3relF = d3*er3;
rG3relO = rErelO + rFrelE + rG3relF;

vG1relO = d1*theta1dot*etheta1;     %first link
%second link
vErelO = L1*theta1dot*etheta1;
vG2relE = d2*theta2dot*etheta2;
vG2relO = vErelO + vG2relE;  
%third link
vFrelE = L2*theta2dot*etheta2;
vG3relF = d3*theta3dot*etheta3;
vG3relO = vErelO + vFrelE + vG3relF;


%Kinetic Energy
Ek1 = 1/2*m1*dot(vG1relO, vG1relO) + 1/2*I1*theta1dot^2;
Ek2 = 1/2*m2*dot(vG2relO, vG2relO) + 1/2*I2*theta2dot^2;
Ek3 = 1/2*m3*dot(vG3relO, vG3relO) + 1/2*I3*theta3dot^2;
Ek = Ek1 + Ek2 + Ek3;

%Potential Energy
Ep1 = -m1*g*rG1relO;    %Negative because 0 level at origin O
Ep2 = -m2*g*rG2relO;
Ep3 = -m3*g*rG3relO;
Ep = Ep1(1) + Ep2(1) + Ep3(1);

%Setting up the Lagrange Equation
L  = Ek - Ep;

thetas = [theta1, theta2, theta3];
thetaDots = [theta1dot, theta2dot, theta3dot];
thetaDoubleDots = [theta1doubledot, theta2doubledot, theta3doubledot];

%d/dt(dL/dxDot) - dL/dx = 0

dLdthetaDots = jacobian(L, thetaDots);
%the first three terms are the left term of the equation (chain rule)
eqns =   jacobian(dLdthetaDots,thetas)*thetaDots'  ...
       + jacobian(dLdthetaDots,thetaDots)*thetaDoubleDots' ...
       + jacobian(dLdthetaDots,t)               ...
       - jacobian(L,thetas)';

[M,b] = equationsToMatrix(eqns,thetaDoubleDots); 
M_L = simplify(M);
b_L = simplify(b);

matlabFunction(M_L, 'file', 'rhsStuffmassMatrixLagrange');
matlabFunction(b_L, 'file', 'rhsstuffbVectorLagrange');