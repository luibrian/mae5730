%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%Parameters
p.gamma = 20*pi/180;            %in rads/s
p.m = 1;                        %mass of rolling cylinder
p.R = 3;                        %radius of cylinder
p.d = 0.5;                      %offset of CM from center
p.I = 1;                        %moment of inertia
p.g = 1;                        %gravity

%timespan
dur = 10;
npoints = 201;
tspan = linspace(0, dur, npoints);

theta0 = 0;
thetaDot0 = 1;
s0 = p.R*theta0;
sDot0 = p.R*thetaDot0;
z0 = [theta0; thetaDot0; s0; sDot0]; %s sDot theta thetaDot

%Doing ODE45
f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);

thetaArray = zArray(:,1);
sArray = zArray(:,3);

figure(1);
plot(tArray, thetaArray);
figure(2);
plot(tArray, sArray);


function zdot = rhs(z,p) 
    %unpacks all the parameters into rhs function using the AMB method
    %names = fieldnames(p);
    %for i = 1:length(names)
    %    eval([names{i} '= p.' names{i} ';']);
    %end
    
    gamma = p.gamma;
    m = p.m;
    R = p.R;
    d = p.d;
    I = p.I;
    g = p.g;
    
    theta = z(1);
    thetaDot = z(2);
    s = z(3);
    sDot =z(4);
    
    b = rhsStuffbVector(R,d,g,gamma,m,thetaDot,theta);
    A = rhsStuffMassStuff(I,R,d,m,theta);
    q =  A\b;
    
    sDoubleDot = q(1);
    thetaDoubleDot = q(2);
    
    zdot = [thetaDot; thetaDoubleDot; sDot; sDoubleDot];
end
