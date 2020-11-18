%Brian Lui
%bl569
%MAE 5730
%HW due 9/12/18
%Cannon ball2- #11

%clear all;
close all;

%specified constants
p.theta = pi/4; theta = p.theta;    %angle
p.v0 = 1; v0 = p.v0;                %[m/s]   Initial velocity
p.g = 1; g = p.g;                   %[m/s^2] Gravity
p.m = 1; m = p.m;                   %[kg]    Mass
p.c = 1; c = p.c;                   %[kg/s]  Drag constant

dur = 20;             %[s]
npoints = 10001;       %number of points
tspan = linspace(0, dur, npoints);

%Initial conditions
%[x(0); xDot(0); y(0); yDot(0); drag work]
%starts with initial velocity, and drag work is 0 because no change in
%position yet
vx0 = p.v0*cos(p.theta);
vy0 = p.v0*sin(p.theta);
z0 = [0; vx0; 0; vy0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan,z0);
xArray = zArray(:,1);
yArray = zArray(:,3);

plot(xArray,yArray);
title('Cannon ball (theta=pi/4, v0=1 m/s, g=1 m/s^2, m=1 kg, c=1kg/s')
xlabel('X position (m)')
ylabel('Y Position (m)')
xlim([0 1])
ylim([-.1 .2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part b. Comparing work done by drag force and change in total energy.

%changing the step sizes
npointsArray = [101 1001 5001 10001 50001 100001];
KEminusWork = zeros(length(npointsArray),1);

for i = 1:length(npointsArray)
    %doing ode45
    tspan = linspace(0, dur, npointsArray(i));
    relTol = 1e-10;
    absTol = 1e-10;
    options = odeset('RelTol', relTol,'AbsTol', absTol);
    
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan,z0);
    
    xArray = zArray(:,1);
    vxArray = zArray(:,2);
    yArray = zArray(:,3);
    vyArray = zArray(:,4);
    work = zArray(:,5);
    
    %finding the index when y=0
    belowGround = find(yArray<0);
    indexOnGround = belowGround(1) -1;

    %Finding the change in KE. It is negative because the work done by drag
    %force is negative
    speedInitial = sqrt(vx0^2 + vy0^2);
    speedFinal = sqrt(vxArray(indexOnGround).^2 + vyArray(indexOnGround).^2);
    deltaKEtot = 1/2*(p.m)*(speedFinal^2 - speedInitial^2);
    
    KEminusWork(i) = deltaKEtot - work(indexOnGround);
end


%There's a steady state difference between drag work and delta KE because
%velocity isn't taken exactly when the ball is back on the ground
figure(2)
plot(npointsArray, abs(KEminusWork))
title('\Delta KE - Work_{drag}')
xlabel('Number of points')
ylabel('Difference between \Delta KE and Drag Work')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function zDot = rhs(z, p)
    %cannon ball dynamics
    x =z(1); xDot = z(2); y = z(3); yDot = z(4); work = z(5);
    g = p.g; m = p.m; c = p.c;
    
    speed = sqrt(xDot^2 + yDot^2);
    
    xDoubleDot = -c/m * speed * xDot;
    yDoubleDot = -c/m*speed*yDot -g;
    workDot = -c*speed*(x*xDot + y*yDot);
    
    zDot = [xDot; xDoubleDot; yDot; yDoubleDot; workDot];
end
