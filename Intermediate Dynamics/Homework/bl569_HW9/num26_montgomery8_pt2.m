%Brian Lui
%MAE 5730
%HW9
%10/27/18

clear all
close all

%26. Want to find a periodic solution by optimizing both the initial
% conditions and simulation time so that the error between initial and
% final state goes to 0. Simulate for 5+ periods to demonstrate that your
% solution is periodic.

%Reasonable guess of initial conditions of a periodic solution. Initial
%guess of period is 8 (time units).
%Initial Conditions
x1_0 = -0.755;
y1_0 = 0.355;
x2_0 = 1.155;
y2_0 = -0.0755;
x3_0 = -0.4055;
y3_0 = -0.3055;

vx1_0 = 0.9955;
vy1_0 = 0.07855;
vx2_0 = 0.1055;
vy2_0 = 0.4755;
vx3_0 = -1.1055;
vy3_0 = -0.5355;
dur0 = 8;

z0 = [x1_0;y1_0;x2_0;y2_0;x3_0;y3_0; vx1_0;vy1_0;vx2_0;vy2_0;vx3_0;vy3_0];

%parameters
p.m1 = 1; p.m2 = 1; p.m3 = 1; p.G = 1;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
g.npoints = 201;         %how many points to put on plot.

%Plot of initial before it has been optimized
tspan = linspace(0,dur0, g.npoints);

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0, options);

x1Array = zArray(:,1);
y1Array = zArray(:,2);
x2Array = zArray(:,3);
y2Array = zArray(:,4);
x3Array = zArray(:,5);
y3Array = zArray(:,6);

figure(1);
plot(x1Array,y1Array, 'b');
hold on;
plot(x2Array,y2Array, 'g');
plot(x3Array,y3Array, 'r--');


%Now let's see if we can get one after fmincon works magic
z0_mincon = [z0; dur0];

h = @(q) toMin(q,p,g);
%Want the period to be at least 15
A = [];
b = [];   
Aeq = [0 0 0 0 1 0 0 0 0 0 0 0 0;...    5x3=0
    0 0 0 0 0 1 0 0 0 0 0 0 0;...       %y3 = 0
    1 0 1 0 0 0 0 0 0 0 0 0 0;...       %x2=-x1
    0 1 0 1 0 0 0 0 0 0 0 0 0];         %y2 = -y1
beq = [0 0 0 0];                    %x3=0, y3=0
lb = [-3 -3 -3 -3 -3 -3,...     %x1,y1,x2,y2,x3,y3
    0, 0, -inf, -inf, -inf, -inf, 7];    %lower bound
ub =[0, 3, 3, 3, 3, 3,...
    inf, inf, inf, inf, inf, inf, inf];     %upper bound
nonlcon = [];
optionsfmincon = optimoptions('fmincon', 'ConstraintTolerance', 1e-14,...
    'StepTolerance', 1e-14);
[bestIC, fval] = fmincon(h, z0_mincon, A, b, Aeq, beq, lb, ub, nonlcon, optionsfmincon);


%Plotting it now with the best IC
tspanBest1Per = linspace(0, bestIC(13), g.npoints);
[tArrayBest, zArrayBest] = ode45(f, tspanBest1Per, bestIC(1:12));
x1ArrayBest = zArrayBest(:,1);
y1ArrayBest = zArrayBest(:,2);
x2ArrayBest = zArrayBest(:,3);
y2ArrayBest = zArrayBest(:,4);
x3ArrayBest = zArrayBest(:,5);
y3ArrayBest = zArrayBest(:,6);

figure(2);
plot(x1ArrayBest, y1ArrayBest);
hold on;
plot(x2ArrayBest, y2ArrayBest);
plot(x3ArrayBest, y3ArrayBest);

%Plotting it now with the best IC, 5 periods
tspanBest5Per = linspace(0, 5*bestIC(13), g.npoints);
[tArrayBest, zArrayBest5] = ode45(f, tspanBest5Per, bestIC(1:12));
x1ArrayBest5 = zArrayBest5(:,1);
y1ArrayBest5 = zArrayBest5(:,2);
x2ArrayBest5 = zArrayBest5(:,3);
y2ArrayBest5 = zArrayBest5(:,4);
x3ArrayBest5 = zArrayBest5(:,5);
y3ArrayBest5 = zArrayBest5(:,6);

figure(3);
plot(x1ArrayBest5, y1ArrayBest5);
hold on;
plot(x2ArrayBest5, y2ArrayBest5);
plot(x3ArrayBest5, y3ArrayBest5);




function zDot = rhs(z,p)
    %This contains the equations of motion for the 3 mass system.
    m1 = p.m1; m2 = p.m2; m3 = p.m3; G = p.G;
    
    x1 = z(1); y1 = z(2); x2 = z(3); y2 = z(4); x3 = z(5); y3 = z(6);
    x1Dot = z(7); y1Dot = z(8); 
    x2Dot = z(9); y2Dot = z(10); 
    x3Dot = z(11); y3Dot = z(12);
    
    el_12 = sqrt((x2-x1)^2 + (y2-y1)^2);
    el_13 = sqrt((x3-x1)^2 + (y3-y1)^2);
    el_23 = sqrt((x3-x2)^2 + (y3-y2)^2);
    
    %equations of motions. Only forces are the gravitational force between
    %them
    x1DoubleDot = G*m2/el_12^3 *(x2-x1) + G*m3/el_13^3 * (x3-x1);
    y1DoubleDot = G*m2/el_12^3 *(y2-y1) + G*m3/el_13^3 * (y3-y1);
    x2DoubleDot = G*m1/el_12^3 *(x1-x2) + G*m3/el_23^3 * (x3-x2);
    y2DoubleDot = G*m1/el_12^3 *(y1-y2) + G*m3/el_23^3 * (y3-y2);
    x3DoubleDot = G*m1/el_13^3 *(x1-x3) + G*m2/el_23^3 * (x2-x3);
    y3DoubleDot = G*m1/el_13^3 *(y1-y3) + G*m2/el_23^3 * (y2-y3);
    
    zDot = [x1Dot; y1Dot; x2Dot; y2Dot; x3Dot; y3Dot;...
        x1DoubleDot; y1DoubleDot;...
        x2DoubleDot; y2DoubleDot;...
        x3DoubleDot; y3DoubleDot];
end



function minFunc = toMin(q,p,g)
    %This contains the equation that needs to be minimized
    %Want to find a periodic solution by optimizing both the initial
    %conditions and simulation time so that the error between initial and
    %final state goes to 0.
    
    %This function gets put in fmincon
    
    %Needs to do ode45 to get the position and velocity of points after
    %simulation of time dur. 
    
    %q contains the state of the masses (x1,y1,x2,y2,x3,y3,
    %vx1,vy1,vx2,vy2,vx3,vy3) and time for period. It's massive because 3
    %points
    x1  = q(1);  y1  = q(2); 
    x2  = q(3);  y2  = q(4);
    x3  = q(5);  y3  = q(6);
    vx1 = q(7);  vy1 = q(8);
    vx2 = q(9);  vy2 = q(10);
    vx3 = q(11); vy3 = q(12);
    period = q(13);
    
    %stuff required for the ode45 function 
    z0_It = [x1 ;y1;...
             x2 ;y2;...
             x3 ;y3;...
             vx1;vy1;...
             vx2;vy2;...
             vx3;vy3];       %Initial condition
         
    npoints = g.npoints;

    tspan = linspace(0, period, npoints);
    
    
    %options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan, z0_It);
    
    %Extracting the vectors out after ode45
    x1Array     = zArray(:,1);
    y1Array     = zArray(:,2);
    x2Array     = zArray(:,3);
    y2Array     = zArray(:,4);
    x3Array     = zArray(:,5);
    y3Array     = zArray(:,6);
    vx1Array    = zArray(:,7);
    vy1Array    = zArray(:,8);
    vx2Array    = zArray(:,9);
    vy2Array    = zArray(:,10);
    vx3Array    = zArray(:,11);
    vy3Array    = zArray(:,12);
    
    %Error between initial position and final position
    x1Error = abs(x1Array(1) - x1Array(end));
    y1Error = abs(y1Array(1) - y1Array(end));
    x2Error = abs(x2Array(1) - x2Array(end));
    y2Error = abs(y2Array(1) - y2Array(end));
    x3Error = abs(x3Array(1) - x3Array(end));
    y3Error = abs(y3Array(1) - y3Array(end));
    
    vx1Error = abs(vx1Array(1) - vx1Array(end));
    vy1Error = abs(vy1Array(1) - vy1Array(end));
    vx2Error = abs(vx2Array(1) - vx2Array(end));
    vy2Error = abs(vy2Array(1) - vy2Array(end));
    vx3Error = abs(vx3Array(1) - vx3Array(end));
    vy3Error = abs(vy3Array(1) - vy3Array(end));
    

    
    %Condition to satisfy for periodic function
    %In theory: 
    %x1Error = x2Error = x3Error
    %y1Error = y2Error = y3Error
    %vx1Error = vx2Error = vx3Error
    %vy1Error = vy2Error = vy3Error
    %but there may be small differences due to numerical errors?
    minPos = 1000*(x1Error^2 + y1Error^2 +...
        x2Error^2 + y2Error^2 + x3Error^2 + y3Error^2);
    %minPos = x1Error + y1Error + x2Error + y2Error + x3Error + y3Error;
    minVel = 100*(vx1Error + vy1Error + vx2Error + vy2Error + vx3Error + vy3Error);
    minFunc = minPos + minVel
end

