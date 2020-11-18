%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%Number 29. Bead on a parabolic wire. Frictionless point mass bead sliding
%on a rigid wire on the curve y= cx^2 with gravity in the -y direction. 

p.m = 1; 
p.g = 1;
p.c = 1;

%setting up tspan
dur = 20;
npoints = 200;
tspan = linspace(0, dur, npoints);

%initial condition
z0_NE = [2; 0]; 
z0_DAE = [2; 0; p.c*2^2; 0];
z0_Lagrange = [2; 0];

%Setting tolerance
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%Doing the ode45
f_NE = @(t,z) rhsNE(z,p);
f_DAE = @(t,z) rhsDAE(z,p);
f_Lagrange = @(t,z) rhsLagrange(z,p);

[tArrayNE, zArrayNE] = ode45(f_NE, tspan, z0_NE);
[tArrayDAE, zArrayDAE] = ode45(f_DAE, tspan, z0_DAE);
[tArrayLagrange, zArrayLagrange] = ode45(f_Lagrange, tspan, z0_Lagrange);

%Extracting the stuff out of the array
xArrayNE = zArrayNE(:,1);
yArrayNE = p.c.* xArrayNE.^2;

xArrayDAE = zArrayDAE(:,1);
yArrayDAE = zArrayDAE(:,3);
%[xArrayDAE, zArrayDAE(:,2), yArrayDAE]

xArrayLagrange = zArrayLagrange(:,1);
yArrayLagrange = p.c.* xArrayLagrange.^2;

%Plotting the stuff
figure(1);
plot(xArrayNE, yArrayNE);
figure(2);
plot(xArrayDAE, yArrayDAE);
figure(3);
plot(xArrayLagrange, yArrayLagrange);


%Checking errors
figure(4);
subplot(2,1,1);
plot(tArrayNE, xArrayNE-xArrayLagrange);
xlabel('time');
ylabel('x Position(NE-Lagrange)');
title('x Position(NE-Lagrange) vs time');
subplot(2,1,2);
plot(tArrayNE, xArrayNE-xArrayDAE);
title('x Position(NE-DAE) vs time');
xlabel('time');
ylabel('x Position(NE-DAE) ');

figure(5);

%Animation of one of them
figure(6);
for i = 1:length(tspan)
   plot(xArrayDAE(i), yArrayDAE(i), '.');
   axis([-5 5 0 10]);
   pause(0.0001);
end




function zdot = rhsNE(z,p)
    %EOM derived through Newton Euler
    g = p.g; c = p.c;      %parameters
    x = z(1); xdot = z(2);          %state
    
    xdoubledot_NE = rhs_NE(c,g,x,xdot);
    %xdoubledot_NE = -2*c*x *(g+2*c*xdot^2)/(1+4*c^2*x^2);
    zdot = [xdot; xdoubledot_NE];
end


function zdot = rhsDAE(z,p)
    %EOM derived through DAE
    m = p.m; g = p.g; c = p.c;                      %parameters
    x = z(1); xdot = z(2); y = z(3); ydot = z(4);   %state
    
    q_DAE = rhs_Matrix_DAE(c,g,m,x,xdot,y);
    xdoubledot = q_DAE(1);
    ydoubledot = q_DAE(2);
    
    zdot = [xdot; xdoubledot; ydot; ydoubledot];
end

function zdot = rhsLagrange(z,p)
    %EOM derived through Lagrange
    g = p.g; c = p.c;      %parameters
    x = z(1); xdot = z(2);          %state
    
    xdoubledot_Lagrange = rhs_Lagrange(c,g,x,xdot);
    
    zdot = [xdot; xdoubledot_Lagrange];
end