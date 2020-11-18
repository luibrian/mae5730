%Brian Lui
%MAE 5730- Intermediate Dynamics
%HW due 8/30/18

clear all;
close all;

%1. ODE & animation practice. 
%Take a simple set of ODEs. Use a set you like,e.g.,
%harmonic oscillator, non-linear pendulum, the Lorentz system (look it up on the
%internet). Solve these numerically 3 ways (see below), and understand the accuracy.
%The goal is that, by the time you hand in the homework, you can write and debug
%the assignment on your own without looking up anything (outside of trivial syntax
%things).

%I'm going to use a non-linear pendulum
%Equation of motion: thetaDoubleDot + (g/L)sin(theta) = 0

%Part a.
%Method 1: as simply as possible, without ODE45, calling functions or anything
%like that. A single function or script with no function calls (ok, plotting calls
%are ok).

%Constants
g = 9.8;        %[m/s^2]    Gravity
L = 3;          %[m]        Length of pendulum
dur = 20;       %[s]
npoints = 10001;
dt = dur/npoints;   %Timestep size
tspan = linspace(0,dur,npoints);

%initializing the variables of interest
time = zeros(npoints+1,1);
theta = zeros(npoints+1,1);
thetaDot = zeros(npoints+1,1);
thetaDoubleDot = zeros(npoints+1,1);

%initial conditions
theta(1) = 0; 
thetaDot(1) = 0.1;

%getting the values of the theta after every time step
for i = 1:npoints
    thetaDoubleDot(i) = -g/L *sin(theta(i));
    time(i+1) = time(i) + dt;
    theta(i+1) = theta(i) + thetaDot(i)*dt;
    thetaDot(i+1) = thetaDot(i) + thetaDoubleDot(i)*dt;
end

figure(1)
plot(time, theta)



%Part b.
%With your own Euler solver function. Your main program should call your Euler
%solver. Your Euler solver should call a RHS (Right Hand Side) function.

%Uses constants from part a. Puts it in parameters of structure p. 
p.g = g;
p.L = L;
timeSpan.dur = 20;
timeSpan.npoints = 1001;

%initial conditions
theta0 = 0; 
thetaDot0 = 0.1;
z0 = [theta0; thetaDot0];

%does both part b and d
figure(2);
nPointsArray = [2, 11, 101, 1001, 10001, 100001];
for n = 1:length(nPointsArray)
    timeSpan.npoints = nPointsArray(n);
    [tArray, zArray] = eulermethod(timeSpan,z0,p);
    thetaArray = zArray(:,1);
    subplot(2,3,n);
    plot(tArray, thetaArray);
    title(['Subplot ', num2str(n), ': ', num2str(nPointsArray(n)), ' points']);
end



%Part c. With ODE45
%Uses the parameters from part a/b, the rhs equation from b
figure(3);
tolerances = [1e-10 1e-8 1e-6 1e-4 1e-2 1];
for i = 1:length(tolerances)
    relTol = tolerances(i);
    absTol = tolerances(i);
    options = odeset('RelTol', relTol,'AbsTol', absTol);

    f = @(t,z)  rhs(z,p); 
    [tArrayODE45, zArrayODE45] = ode45(f, tspan, z0, options);
    thetaArrayODE45 = zArrayODE45(:,1);
    
    subplot(2,3,i);
    plot(tArrayODE45, thetaArrayODE45);
    titleString = ['Relative and Absolute tolerance of ', num2str(tolerances(i))];
    title(titleString);
end


%Part d.
%Using (b) solve the equations many times with progressively smaller step size,
%down to the smallest size you have patience for, and up to the largest size that
%isn't crazy. As sensibly as possible, compare the results and use that comparison
%to estimate the accuracy of each solution. You should be able to find a method
%to estimate the accuracy of a numerical solution even without knowing the exact
%solution.

%I looped through various step sizes and plotted them in part b. The smaller
%the step size, the more accurate the graph was. If you increase the number
%of points but the graph of the solution does not noticeably change, then
%you can say with a relatively high degree of certainty that the numerical
%solution is correct.



%Part e.
%Using ODE45, solve the equations with various accuracies (use 'reltol' and 'ab-
%stol'). Does Matlab do a good job of estimating its own accuracy? Use suitable
%plots to make your point.

%Matlab does a pretty good job estimating its own accuracy. Even with large
%tolerances, it does a good job plotting the harmonic function.




function [tArray, zArray] = eulermethod(timeSpan, z0, p)
%for this function, need to pass in information about:
%time to create tArray
%initial theta and thetaDot condtion within z0
%p for other relevant parameters
%this function is used in part (b)

    dur = timeSpan.dur;
    npoints = timeSpan.npoints;
    tArray = linspace(0, dur, npoints);
    dt = dur/npoints;
    
    zArray = zeros(npoints, 2);
    zArray(1,:) = z0;
    
    for i = 1:(npoints-1)
        zArray(i+1,:) = zArray(i,:) + dt*rhs(zArray(i,:),p)';
        %the next element is the current element * change in time
    end
end


function zdot = rhs(z,p)
   g = p.g;
   L = p.L;
   theta = z(1);
   thetaDot = z(2);
   
   thetaDoubleDot = -g/L * sin(theta);
   zdot = [thetaDot; thetaDoubleDot];
end
