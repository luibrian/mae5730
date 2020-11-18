%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%a. 
%Draw a box that is long in the x direction and narrow in 
%the z direction so the moment of inertia in the reference
%configuration could plausibly be:
%I = [1 0 0; 0 2 0; 0 0 3];

clear all;
close all;

%Setting up I to have the desired stuff
m = 1;
b = sqrt(24/m) + 0.5;          %length of box base
w = sqrt(12/m) + 0.5;          %length of box width
h = 0.5;                       %length of box height

%Moment of inertia matrix
I0 = [w^2+h^2,       0,        0;...
           0, b^2+h^2,        0;...
           0,       0, b^2+w^2] * m/12;
I0 = round(I0);

%input parameters
p.I0 = I0;      %initial inertia of object
p.I1 = I0(1,1); p.I2 = I0(2,2); p.I3 = I0(3,3);
p.M = 0;        %net moment acting on rigid object


opt = odeset('Reltol', 1e-10, 'Abstol', 1e-10);

%b.
%Set up numerical solutions of the general non-linear 
%free motion of such an object. That is, given initial 
%values for the body-fixed angular velocity components, 
%find the body-fixed angular velocities as a function 
%of time.- This is in the rhs function.

%c,d.
%As an initial condition, consider a small deviation 
%to rotation about the x' axis, y' axis, z' axis. 

dur = 10;
npoints = 501;
tspan = linspace(0, dur, npoints);

dev_initial = [1e-2 1e-2 1e-2]';
%small deviations to rotation about x',y',z' axis
x_IC = [1 0 0]' + dev_initial;  %case 1
y_IC = [0 1 0]' + dev_initial;  %case 2
z_IC = [0 0 1]' + dev_initial;  %case 3

f = @(t,z) rhs(z,p);

[tArray_x, zArray_x] = ode45(f, tspan, x_IC, opt);
[tArray_y, zArray_y] = ode45(f, tspan, y_IC, opt);
[tArray_z, zArray_z] = ode45(f, tspan, z_IC, opt);

%Plotting the small deviations of omega over time
%case 1
w1_smallDev_x_Array = zArray_x(:,1) - ones(length(tArray_x));
w2_smallDev_x_Array = zArray_x(:,2) - zeros(length(tArray_x));
w3_smallDev_x_Array = zArray_x(:,3) - zeros(length(tArray_x));

%case 2
w1_smallDev_y_Array = zArray_y(:,1) - zeros(length(tArray_y));
w2_smallDev_y_Array = zArray_y(:,2) - ones(length(tArray_y));
w3_smallDev_y_Array = zArray_y(:,3) - zeros(length(tArray_y));

%case 3
w1_smallDev_z_Array = zArray_z(:,1) - zeros(length(tArray_z));
w2_smallDev_z_Array = zArray_z(:,2) - zeros(length(tArray_z));
w3_smallDev_z_Array = zArray_z(:,3) - ones(length(tArray_z));


%Plotting case 1
figure(1);
subplot(3,1,1);
plot(tArray_x, w1_smallDev_x_Array);
title('w1 for small deviations in rotation about x axis')
subplot(3,1,2);
plot(tArray_x, w2_smallDev_x_Array);
title('w2 for small deviations in rotation about x axis')
subplot(3,1,3);
plot(tArray_x, w3_smallDev_x_Array);
title('w3 for small deviations in rotation about x axis')


%Plotting case 2
figure(2);
subplot(3,1,1);
plot(tArray_y, w1_smallDev_y_Array);
title('w1 for small deviations in rotation about y axis')
subplot(3,1,2);
plot(tArray_y, w2_smallDev_y_Array);
title('w2 for small deviations in rotation about y axis')
subplot(3,1,3);
plot(tArray_y, w3_smallDev_y_Array);
title('w3 for small deviations in rotation about y axis')


%Plotting case 3
figure(3);
subplot(3,1,1);
plot(tArray_z, w1_smallDev_z_Array);
title('w1 for small deviations in rotation about z axis')
subplot(3,1,2);
plot(tArray_z, w2_smallDev_z_Array);
title('w2 for small deviations in rotation about z axis')
subplot(3,1,3);
plot(tArray_z, w3_smallDev_z_Array);
title('w3 for small deviations in rotation about z axis')


%f.
%Using some non-trivial initial condition, show that 
%the solution to these equations agrees with the 
%solutions from the full equations with rotations. 
%How to compare? Use the fixed frame w components and 
%R as a function of time, from the previous problem, 
%to find the body-fixed components of w as a function 
%of time.

%ic = [1 1 1]' + dev_initial;
ic = [0.01, 0.1, 0]';

%Finding out the motion by same method above
[tArray_ic_b, zArray_ic_b] = ode45(f, tspan, ic);
w1_icb_Array = zArray_ic_b(:,1);
w2_icb_Array = zArray_ic_b(:,2);
w3_icb_Array = zArray_ic_b(:,3);

%Finding out the motion by num47 
R_0 = eye(3);                      %Initally not rotating
%State: R_0 as 9x1, ic as 3x1
z0 = [reshape(R_0, [9,1]); ic];

g = @(t,z) rhs47(z,p);
[tArray_ic_f, zArray_ic_f] = ode45(g, tspan, z0, opt);
w1_icf_Array = zArray_ic_f(:,10);
w2_icf_Array = zArray_ic_f(:,11);
w3_icf_Array = zArray_ic_f(:,12);

figure(4);
subplot(3,1,1);
plot(tArray_ic_b, w1_icf_Array- w1_icb_Array)
title('w1 for w0 = [0.01, 0.1, 0]');
subplot(3,1,2);
plot(tArray_ic_b, w2_icf_Array- w2_icb_Array)
title('w2 for w0 = [0.01, 0.1, 0]');
subplot(3,1,3);
plot(tArray_ic_b, w3_icf_Array- w3_icb_Array)
title('w3 for w0 = [0.01, 0.1, 0]');
%Result: It's very close together, with the numerical error accumulating
%over time







function zdot = rhs(z,p)
    %Done with Euler equations
    w1 = z(1); w2 = z(2); w3 = z(3);
    I1 = p.I1; I2 = p.I2; I3 = p.I3;
    
    w1dot = (I2-I3)/I1 *w2 *w3;
    w2dot = (I3-I1)/I1 *w1 *w3;
    w3dot = (I1-I2)/I1 *w1 *w2;
    
    zdot = [w1dot w2dot w3dot]';
end


function zdot = rhs47(z, p)
    %Doing rhs using question 47. 
    %Does doing the full I, omega tensors
    I0 = p.I0;
    M = p.M;
    %extracting the state
    R = reshape(z(1:9,1), [3,3]);
    omegaVector = z(10:12,1);
    
    %This is the skew symmetric version of omegaVector
    omegaTensor = [0, -omegaVector(3),  omegaVector(2);...
      omegaVector(3),               0, -omegaVector(1);...
     -omegaVector(2),  omegaVector(1),              0];
    
    I = R * I0 * R';
    omegaVectorDot = inv(I) * (eye(3,1)*M - cross(omegaVector, I*omegaVector));
    Rdot = omegaTensor * R;
    zdot = [reshape(Rdot, [9,1]); omegaVectorDot];
end
