%Brian Lui
%bl569
%MAE5730

%Deriving the EOM for the eccentric cylinder 

syms s theta sDot thetaDot sDoubleDot thetaDoubleDot real
syms Ff Fn real %friction force, normal force
syms m g I real %mass, gravity, moment of inertia
syms R d gamma real % radius, distance off center, slope

%fixed frame (ex ey ez), fixed frame aligned with slope (elamb, en), 
%rotating frame (er etheta)
syms ex ey ez elamb en er etheta real
%vectors of interest- O is inertial frame point, G is CM, P is point in
%contact with surface
syms rGrelO vGrelO aGrelO real %combined vectors
syms rGrelP vGrelP aGrelP real %CM relative to point on ground
syms rPrelO vPrelO aPrelO real %point on ground relative to inertial

%Defining the frames
%Inertial frame
ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

%Inertial frame aligned with rotating cylinder
elamb = cos(gamma)*ex - sin(gamma)*ey;
en = cross(ez, elamb);

%Rotating frame, er always pointing from center C to CM G
er = cos(theta)*elamb - sin(theta)*en;
etheta = cross(ez,er);

%Important pos/vel/accel vectors
rPrelO = s*elamb;
rGrelP = R*en + d*er;
rGrelO = rGrelP + rPrelO;

vPrelO = sDot*elamb; 
vGrelP = - d*thetaDot*etheta;
vGrelO = vGrelP + vPrelO;

aPrelO = sDoubleDot*elamb;
aGrelP = - d*thetaDot^2*er - d*thetaDoubleDot*etheta;
aGrelO = aPrelO + aGrelP;



%Forces on system
gravForce = -m*g*ey;
fricForce = -Ff*elamb;
normForce = Fn*en;


%Linear momentum balance
lmb = m*aGrelO - fricForce - normForce - gravForce;
eqn1 = dot(lmb, elamb);
eqn2 = dot(lmb, en);

%No slip equation
eqn3 = sDoubleDot - R*thetaDoubleDot;

%Angular momentum balance
moments = dot(cross(rGrelP, fricForce),ez);
amb_rhs = dot(cross(rGrelP, m*aGrelO) + I*thetaDoubleDot*ez,ez);
eqn4 = amb_rhs - moments;

eqns = [eqn1; eqn2; eqn3; eqn4];
unknowns = [sDoubleDot, thetaDoubleDot, Ff, Fn];

[A,b] = equationsToMatrix(eqns, unknowns);
A = simplify(A);
b = simplify(b);

matlabFunction(A, 'file', 'rhsStuffMassStuff');
matlabFunction(b, 'file', 'rhsStuffbVector');
