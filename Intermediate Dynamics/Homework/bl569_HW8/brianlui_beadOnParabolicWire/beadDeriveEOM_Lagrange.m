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

dist = sqrt(x^2 + y^2);
cosineT = x/dist;
sineT = y/dist; 

et = cosineT*e1 + sineT * e2;
en = cross(e3,et);

%Constraint- there is only 1 equation;
y = c*x^2;
ydot = 2*c*x*xdot;
ydoubledot = 2*c*xdot^2 + 2*c*xdoubledot*x;

%Position vector
rPrelO = x*e1 + y*e2;
vPrelO = xdot*e1 + ydot*e2;
aPrelO = xdoubledot*e1 + ydoubledot*e2;


%Doing it using Lagrange
Ek = 1/2*m*dot(vPrelO,vPrelO);
Ep = m*g*y;

L = Ek - Ep;
dLdx = jacobian(L,x);
dLdxdot = jacobian(L,xdot);
%d(dLdxdot)/dt
ddLdxdotdt = simplify(jacobian(dLdxdot,xdot)*xdoubledot) +...
    jacobian(dLdxdot,x)*xdot + ...
    jacobian(dLdxdot,t);


eqn_Lagrange = simplify(ddLdxdotdt - dLdx);

xdoubledot_Lagrange = simplify(solve(eqn_Lagrange, xdoubledot));

matlabFunction(xdoubledot_Lagrange,'file','rhs_Lagrange');
