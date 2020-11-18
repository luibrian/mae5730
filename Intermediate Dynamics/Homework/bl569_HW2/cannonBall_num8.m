%Brian Lui
%bl569
%MAE 5730
%HW due 9/12/18
%Cannon ball- #8

%clear all;
close all;

%specified constants
p.theta = pi/4; theta = p.theta;    %angle
p.v0 = 1; v0 = p.v0;                %[m/s]   Initial velocity
p.g = 1; g = p.g;                   %[m/s^2] Gravity
p.m = 1; m = p.m;                   %[kg]    Mass
p.c = 1; c = p.c;                   %[kg/s]  Drag constant
t = 2;   %This is the time at which to compare the num and analyt solutions

dur = 20;             %[s]
npoints = 10001;       %number of points
tspan = linspace(0, dur, npoints);

%Initial conditions
%[x(0); xDot(0); y(0); yDot(0)] 
z0 = [0; p.v0*cos(theta); 0; p.v0*sin(theta)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part b. Find a numerical solution for theta=pi/4, v0=1, g=1, m=1, c=1
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
ylim([0 .2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part c- compare the numerical and analytical solutions. At t=2, how big is
%the error? How does erro depend on specified tolerances or step sizes?
%analytical solution at time t=2
xAnalyt = v0*cos(theta) * (m/c) * (1 - exp(-c/m*t));
yAnalyt = (-g*m*t/c) - (m/c)*(v0*sin(theta) + (g*m/c))*exp(-c/m*t) +...
    m/c* (v0*sin(theta) + g*m/c);

%changing the relative and absolute tolerance
tolerances = [1e-10 1e-8 1e-6 1e-4 1e-2 1];
errors = zeros(length(tolerances),2);
figure(2)                   %plot of part c with the different tolerances

for i = 1:length(tolerances)
    %generates plots of the EOM at different tolerances
    relTol = tolerances(i);
    absTol = tolerances(i);
    options = odeset('RelTol', relTol,'AbsTol', absTol);
    
    %doing ode45
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan,z0,options);
    
    xPlot = zArray(:,1);
    yPlot = zArray(:,3);
    
    subplot(3,2,i)
    plot(xPlot,yPlot)
    title(['Tolerance = ', num2str(tolerances(i))]);
    xlim([0 1])
    ylim([0 .2])
    
    
    %finding the numeric solution at time t=2, comparing it with
    %the analytical, and sticking in error matrix
    index = find(tArray==t);
    xNum = xPlot(index);
    yNum = yPlot(index);
    
    xError = xAnalyt - xNum;
    yError = yAnalyt - yNum;
    errors(i,1) = xError;
    errors(i,2) = yError;
end

errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use larger and larger v0 to see how trajectory changes
%Extracting x and y positions and plotting (for parts b and d)
v0Array = [1 5 10 15 20 25];
figure(3)
for i = 1:length(v0Array)
    p.v0 = v0Array(i);
    z0 = [0; p.v0*cos(theta); 0; p.v0*sin(theta)];
    
    [tArray, zArray] = ode45(f, tspan,z0);
    
    xPlot = zArray(:,1);
    yPlot = zArray(:,3);
    
    plot(xPlot, yPlot)
    hold on
end
legend('v0=1', 'v0=5', 'v0=10', 'v0=15', 'v0=20', 'v0=25');
xlabel('x position (m)')
ylabel('y position (m)')
xlim([0 20])
ylim([0 20])
hold off

%plotting for v0 = infinity
p.v0 = 10000;
z0 = [0; p.v0*cos(theta); 0; p.v0*sin(theta)];
tspanInf = linspace(0,10000,10001);

[tArray, zArray] = ode45(f, tspanInf,z0);
    
xPlot = zArray(:,1);
yPlot = zArray(:,3);

figure(4)
plot(xPlot, yPlot)
xlabel('x position (m)')
ylabel('y position (m)')
xlim([0 10000])
ylim([0 10000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%part e)
%Find the best launch angle thetaStar for maximizing the range for any
%given v0. As v0 approaches inifinity, what angle does thetaStar tend?
v0Array = [1 5 10 15 20 25];
tspanE = linspace(0,30,1001);   %need a larger time interval because range may increase
thetaArray = linspace(0,90,91); %array from 0 to 90 1 degree at a time. NOTE: degrees
maxRangeAtTheta = zeros(length(thetaArray),1); %each value corresponds to the theta w/ same index
thetaStar = zeros(length(v0Array),1);   %max theta corresponding to v0 in v0Array
for i = 1:length(v0Array)
    p.v0 = v0Array(i);
    
    for j = 1:length(thetaArray)
        p.theta = thetaArray(j);
        z0 = [0; p.v0*cosd(p.theta); 0; p.v0*sind(p.theta)];
        
        [tArray, zArray] = ode45(f, tspanE,z0);
        xArray = zArray(:,1);
        yArray = zArray(:,3);

        indexBelow = find(yArray < 0);  %matrix of indices when it's below ground
        %first index-1 is just before it goes underground
        indexGround = indexBelow(1) - 1;
        maxRangeAtTheta(j) = xArray(indexGround);  %max range for this theta
    end
    maxRangeAtTheta;
    furthestRange = max(maxRangeAtTheta);

    indexOfFurthestRange = find(maxRangeAtTheta==furthestRange);
    thetaStar(i) = indexOfFurthestRange;
end

thetaStar
%As velocity increases, thetastar decreases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function zDot = rhs(z, p)
    %cannon ball dynamics
    x =z(1); xDot = z(2); y = z(3); yDot = z(4);
    g = p.g; m = p.m; c = p.c;
    
    xDoubleDot = -c/m * xDot;
    yDoubleDot = -c/m*yDot -g;
    zDot = [xDot; xDoubleDot; yDot; yDoubleDot];
end
