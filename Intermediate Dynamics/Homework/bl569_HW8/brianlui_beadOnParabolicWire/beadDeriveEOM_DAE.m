%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%Derive the equation of motion for a bead on parabolic wire.
%Derives using Newton Euler method

syms e1 e2 e3 et en real                        %unit vectors
syms x xdot xdoubledot y ydot ydoubledot real   %states
syms m c g t real                               %paramaters
syms N real                                     %constraint force
syms rPrelO vPrelO aPrelO                       %relevant vectors

e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

% dist = sqrt(x^2 + y^2);
% cosineT = x/dist;
% sineT = y/dist; 

speed = sqrt(xdot^2 + ydot^2);
cosineT = xdot/speed;
sineT = ydot/speed;

rPrelO = x*e1 + y*e2;
vPrelO = xdot*e1 + ydot*e2;
aPrelO = xdoubledot*e1 + ydoubledot*e2;

%xPos = rPrelO(1);
xVel = vPrelO(1);

%et = vPrelO/sqrt(xdot^2 + ydot^2)
et = 2*c*x*xdot
%et = cosineT*e1 + sineT*e2
en = sign(xVel)*cross(e3,et)

eqnLMB = m*aPrelO - N*en + m*g*e2

%Doing it using DAE
DAEeqn1 = dot(eqnLMB, e1)
DAEeqn2 = dot(eqnLMB, e2)
DAEeqn3 = ydoubledot - (2*c*x)*xdoubledot - 2*c*xdot^2; %constraint

eqns_DAE = [DAEeqn1, DAEeqn2, DAEeqn3];
unknowns_DAE = [xdoubledot, ydoubledot, N];

[A_DAE, b_DAE] = equationsToMatrix(eqns_DAE, unknowns_DAE)
q_DAE = simplify(A_DAE\b_DAE)

matlabFunction(q_DAE, 'file', 'rhs_Matrix_DAE');


