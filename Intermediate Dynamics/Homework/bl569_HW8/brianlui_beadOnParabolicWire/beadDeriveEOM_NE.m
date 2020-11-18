%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%Derive the equation of motion for a bead on parabolic wire.
%Derives using Newton Euler method

syms e1 e2 e3 et en real                        %unit vectors
syms x xdot xdoubledot y ydot ydoubledot real   %states
syms m c g real                                 %paramaters
syms N real                                     %constraint force
syms rPrelO vPrelO aPrelO                       %relevant vectors

e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

%This is the path constraint
y = 2*c*x^2;
ydot = 2*c*x*xdot;
ydoubledot = 2*c*xdot^2 + 2*c*xdoubledot*x;

rPrelO = x*e1 + y*e2;
vPrelO = xdot*e1 + ydot*e2;
aPrelO = xdoubledot*e1 + ydoubledot*e2;

%et = vPrelO/sqrt(vPrelO(1)^2 + vPrelO(2)^2); also okay
et = vPrelO/norm(vPrelO);
en = cross(e3,et);

eqnLMB = m*aPrelO - N*en + m*g*e2;

%Doing it using Newton Euler
%LMB
eqn = simplify(dot(eqnLMB, vPrelO));
xdoubledot_NE = solve(eqn, xdoubledot);

matlabFunction(xdoubledot_NE, 'file', 'rhs_NE');





