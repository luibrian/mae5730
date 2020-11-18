%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%This derives the EOM for the double pendulum. This does it using Lagrange 
%equations.

%Need the real or else Matlab does some funky stuff with complex conjugates
syms d1 d2 L1 m1 I1 m2 I2 g t real          %parameters (Note d one, not d L)
syms theta1 theta2 theta1dot theta2dot theta1doubledot theta2doubledot real %states
syms e1 e2 e3 er1 etheta1 er2 etheta2 real               %unit vectors
syms rG1relO aG1relO real                                %relevant vectors
syms rG2relO aG2relO rG2relE aG2relE rErelO aErelO real  %relevant vectors
syms MrelO MrelE H1relODot H2relODot H real

%defining unit vectors
e1 = [1 0 0]'; 
e2 = [0 1 0]'; 
e3 = cross(e1,e2); 
er1 = cos(theta1)*e1 + sin(theta1)*e2;
etheta1 = cross(e3, er1);
er2 = cos(theta2)*e1 + sin(theta2)*e2;
etheta2 = cross(e3, er2);


%defining position and velocity vectors of interest
rG1relO = d1*er1;                   %first link

%second link
rErelO  = L1*er1;  
rG2relE = d2*er2;
rG2relO = rErelO + rG2relE;

vG1relO = d1*theta1dot*etheta1;     %first link
%second link
vErelO = L1*theta1dot*etheta1;
vG2relE = d2*theta2dot*etheta2;
vG2relO = vErelO + vG2relE;  



%Kinetic Energy
Ek1 = 1/2*m1*dot(vG1relO, vG1relO) + 1/2*I1*theta1dot^2;
Ek2 = 1/2*m2*dot(vG2relO, vG2relO) + 1/2*I2*theta2dot^2;
Ek = Ek1 + Ek2;

%Potential Energy
Ep1 = -m1*g*rG1relO;    %Negative because 0 level at origin O
Ep2 = -m2*g*rG2relO;
Ep = Ep1(1) + Ep2(1);

%Setting up the Lagrange Equation
L  = Ek - Ep;

thetas = [theta1, theta2];
thetaDots = [theta1dot, theta2dot];
thetaDoubleDots = [theta1doubledot, theta2doubledot];

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