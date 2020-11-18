%Brian Lui
%bl569
%MAE 5730
%HW due 9/12/18
%Cannon ball- #8

clear all;
close all;

%specified constants
p.m = 1;                    %[kg]    Mass
p.a = 1;                    

%Before optimizing
%Initial conditions
%[r(0); theta(0); rDot(0); Dot(0);period] 
r0 =1; theta0 = 0; rDot0 = 1; thetaDot0 = 1;
period0 = 25.2;
z0 = [r0; theta0; rDot0; thetaDot0; period0];

z0_ODE = z0(1:4);
%r(0) cannot be 0 initially- get a divide by 0

dur = period0;          %[s]
npoints = 10001;       %number of points
tspan = linspace(0, dur, npoints);

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0_ODE);
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
axis equal;
axis([-2.5 2.5 -2.5 2.5])

h = @(y) toMinimize(y,p); 
%Want the period to be at least 15
A = [];
b = [];
Aeq = [];
beq = [];   
lb = [0.1 0.1 0.1 0.1 5];    %lower bound
ub =[inf inf inf inf inf];     %upper bound
nonlcon = [];
options = optimoptions('fmincon', 'ConstraintTolerance', 1e-12);
[bestIC, fval] = fmincon(h, z0, A, b, Aeq, beq, lb, ub, nonlcon, options);


%Plotting the best condition
bestICForODE = bestIC(1:4);
tspanBest = linspace(0, bestIC(5), 1001);

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f,tspanBest, bestICForODE);
rArray = zArray(:,1);
thetaArray = zArray(:,2);
xArray = rArray.*cos(thetaArray);
yArray = rArray.*sin(thetaArray);

figure(2)
plot(xArray, yArray);
axis equal
axis([-0.7 0.7 -0.7 0.7])


function zDot = rhs(z,p)
    %assuming central force = -1/r^3 * e_r
    r = z(1); theta = z(2);
    rDot = z(3); thetaDot = z(4);
    a = p.a;
    m = p.m;
    
    centralForce = -a*r^(-1/2);
    
    rDoubleDot = 1/m*(r*thetaDot^2 + centralForce);
    thetaDoubleDot = -2*rDot*thetaDot/r;
    
    zDot = [rDot; thetaDot; rDoubleDot; thetaDoubleDot];
end


function yMinned = toMinimize(y, p)
    %gets put into fmincon
    %fmincon should try to minimize this function
    %y is the state, same state as rhs plus 1 more parameter period
    %q is parameters
    
    %y = [r, theta, rDot, thetaDot, period];
    
    %Goal is to minimize r(end)-r(1) + theta(end)-theta(1) +
    %rDot(end)-rDot(1) + thetaDot(end)-thetaDot(1)
    
    %tspan changes every time because it should 
    dur = y(5);             %[s]
    npoints = 1001;       %number of points
    tspan = linspace(0, dur, npoints);
    
    z0_ODE = [y(1); y(2); y(3); y(4)];
    
    %running ode45
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan, z0_ODE);
    
    rArray = zArray(:,1);
    thetaArray = zArray(:,2);
    rDotArray = zArray(:,3);
    thetaDotArray = zArray(:,4);
    
    x = rArray.*thetaArray;
    y = rArray.*thetaArray;
    xDot = -rArray.*sin(thetaArray).*thetaDotArray + rDotArray.*cos(thetaArray);
    yDot = rDotArray.*sin(thetaArray) + rArray.*cos(thetaArray).*thetaDotArray;
    
    %Minifunctions that I want to minimize: want the initial and final 
    %position and velocity of the points to be the same    
    x1 = rArray(1)*cos(thetaArray(1));
    y1 = rArray(1)*sin(thetaArray(1));
    xEnd = rArray(end)*cos(thetaArray(end));
    yEnd = rArray(end)*(sin(thetaArray(end)));
    xDot1 = -rArray(1)*sin(thetaArray(1)).*thetaDotArray(1) + rDotArray(1).*cos(thetaArray(1));
    yDot1 = rDotArray(1).*sin(thetaArray(1)) + rArray(1)*cos(thetaArray(1))*thetaDotArray(1);
    xDotEnd = -rArray(end)*sin(thetaArray(end)).*thetaDotArray(end) + rDotArray(end).*cos(thetaArray(end));
    yDotEnd = rDotArray(end).*sin(thetaArray(end)) + rArray(end)*cos(thetaArray(end))*thetaDotArray(end);
    
    xErr = abs(xEnd - x1);
    yErr = abs(yEnd- y1);
    xDotErr = abs(xDotEnd - xDot1);
    yDotErr = abs(yDotEnd - yDot1);

    yMinned = 10*xErr^2 + 10*yErr^2 + xDotErr + yDotErr

%     r1Err = abs(rArray(end) - rArray(1));
%     thetaErr = abs(thetaArray(end) - thetaArray(1));
%     r1DotErr = abs(rDotArray(end) - rDotArray(1));
%     thetaDotErr = abs(thetaDotArray(end) - thetaDotArray(1));
%     
%     yMinned = 10*r1Err^2 + thetaErr + r1DotErr + thetaDotErr
end