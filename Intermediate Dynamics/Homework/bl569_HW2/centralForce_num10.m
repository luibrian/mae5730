%Brian Lui
%bl569
%MAE 5730
%HW due 9/12/18
%Cannon ball- #8

%clear all;
close all;

%specified constants
p.g = 1; g = p.g;                   %[m/s^2] Gravity
p.m = 1; m = p.m;                   %[kg]    Mass
p.a = 1;

dur = 200;             %[s]
npoints = 10001;       %number of points
tspan = linspace(0, dur, npoints);

%Initial conditions
%[r(0); theta(0); rDot(0); Dot(0)] 
r0 =1; theta0 = 0;
rDot0 = 1; thetaDot0 = 1;
z0 = [1; 1; 1; 1];
%r(0) cannot be 0 initially- get a divide by 0

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);
rArray = zArray(:,1);
thetaArray = zArray(:,2);
rDotArray = zArray(:,3);
thetaDotArray = zArray(:,4);


xArray = rArray.*cos(thetaArray);
yArray = rArray.*sin(thetaArray);

figure(1);
plot(xArray, yArray);
xlabel('x position');
ylabel('y position');


% fun = @periodicFun;
% x0 = [rArray(1); thetaArray(1); rDotArray(1); thetaDotArray(1)];
% fsolve(fun,x0)

function zDot = rhs(z,p)
    %assuming central force = -1/r^3 * e_r
    r = z(1); theta = z(2);
    rDot = z(3); thetaDot = z(4);
    a = p.a;
    m = p.m;
    
    %%%%%%%%%%%%
    %not cool graphs
    %centralForce = -a/r^3  line
    %centralForce = a/r^3   line
    
    %cool graphs- just comment out centralForce line and pop a new one in!
    %centralForce = -a/r
    %centralForce = -a*sqrt(r);
    %centralForce = -a/r^(1/2);
    %centralForce = -a*r/r^(1/2) + 3/r^(1/2);
    
    %%%%%%%%%%%%
    centralForce = -a*r^(-3/2) - a*r^(-5/2);
    rDoubleDot = 1/m*(r*thetaDot^2 + centralForce);
    thetaDoubleDot = -2*rDot*thetaDot/r;
    
    zDot = [rDot; thetaDot; rDoubleDot; thetaDoubleDot];
end

% function f = periodicFun(state)
%     %periodic function will satisfy the four functions f(1)- f(4)
%     
%     rAtT = state(1); thetaAtT = state(2); 
%     rDotAtT = state(3); thetaDotAtT = state(4);
%     
%     f(1) = rAtT - 1;
%     f(2) = thetaAtT - 1;
%     f(3) = rDotAtT - 1;
%     f(4) = thetaDotAtT - 1;
% end