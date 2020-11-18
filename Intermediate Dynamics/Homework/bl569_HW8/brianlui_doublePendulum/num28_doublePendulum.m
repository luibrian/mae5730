%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%28. Double pendulum. Hinges at origin O and elbow E. Both bars are
%unifirom with same length el = 1 (in some consistent unit system), g=1.
%Neglect all friction and assume there are no joint motors. 

%Parameters
p.L1 = 1; p.L2 = 1; p.d1 = 0.5; p.d2 = 0.5; p.m1 = 1; p.m2 = 1; 
%p.I1 = 1/12*p.m1*p.d1^2; p.I2 = 1/12*p.m2*p.d2^2; 
p.g = 1; p.I1 = 0.2; p.I2 = 0.3;

%setting up tspan
dur = 20;
npoints = 200;
tspan = linspace(0, dur, npoints);

%initial condition
z0 = [pi/2; pi; 0; 0];      %upper arm horizontal, forearm vertical up

%Setting tolerance
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the AMB Method
fAMB = @(t,z) rhsAMB(z,p);
[tArrayAMB, zArrayAMB] = ode45(fAMB, tspan, z0, options);

theta1ArrayAMB = zArrayAMB(:,1);
theta2ArrayAMB = zArrayAMB(:,2);

%Getting the points of interest
xEndArrayAMB = p.L1.*sin(theta1ArrayAMB) + p.L2.*sin(theta2ArrayAMB);
yEndArrayAMB = -(p.L1.*cos(theta1ArrayAMB) + p.L2.*cos(theta2ArrayAMB));

xOAMB = zeros(length(tArrayAMB), 1);
yOAMB = zeros(length(tArrayAMB), 1);

xEArrayAMB = p.L1.*sin(theta1ArrayAMB);
yEArrayAMB = -(p.L1.*cos(theta1ArrayAMB));

%The thing to actually plot
xPlotAMB = [xOAMB, xEArrayAMB, xEndArrayAMB]; 
yPlotAMB = [yOAMB, yEArrayAMB, yEndArrayAMB];
 
figure(1)
for i = 1:length(tArrayAMB)
    plot(xPlotAMB(i,:),yPlotAMB(i,:),'LineWidth', 5)  
    axis equal
    axis ([-3 3 -3 3]);
    title('Animation of Double Pendulum-AMB')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end

figure(2);      %just the path of it
plot(xEndArrayAMB, yEndArrayAMB);
title('Position of the end of 2nd link over time-AMB')
xlabel('x Position')
ylabel('y Position')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the DAE Method
fDAE = @(t,z) rhsDAE(z,p);
[tArrayDAE, zArrayDAE] = ode45(fDAE, tspan, z0, options);

theta1ArrayDAE = zArrayDAE(:,1);
theta2ArrayDAE = zArrayDAE(:,2);

%Getting the points of interest
xEndArrayDAE = p.L1.*sin(theta1ArrayDAE) + p.L2.*sin(theta2ArrayDAE);
yEndArrayDAE = -(p.L1.*cos(theta1ArrayDAE) + p.L2.*cos(theta2ArrayDAE));

xODAE = zeros(length(tArrayDAE), 1);
yODAE = zeros(length(tArrayDAE), 1);

xEArrayDAE = p.L1.*sin(theta1ArrayDAE);
yEArrayDAE = -(p.L1.*cos(theta1ArrayDAE));

%The thing to actually plot
xPlotDAE = [xODAE, xEArrayDAE, xEndArrayDAE]; 
yPlotDAE = [yODAE, yEArrayDAE, yEndArrayDAE];
 
figure(3)
for i = 1:length(tArrayDAE)
    plot(xPlotDAE(i,:),yPlotDAE(i,:),'LineWidth', 5)  
    axis equal
    axis ([-3 3 -3 3]);
    title('Animation of Double Pendulum-DAE')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end

figure(4);      %just the path of it
plot(xEndArrayDAE, yEndArrayDAE);
title('Position of the end of 2nd link over time-DAE')
xlabel('x Position')
ylabel('y Position')

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the Lagrange Method
fLagrange = @(t,z) rhsLagrange(z,p);
[tArrayLagrange, zArrayLagrange] = ode45(fLagrange, tspan, z0, options);

theta1ArrayLagrange = zArrayLagrange(:,1);
theta2ArrayLagrange = zArrayLagrange(:,2);

%Getting the points of interest
xEndArrayLagrange = p.L1.*sin(theta1ArrayLagrange) + p.L2.*sin(theta2ArrayLagrange);
yEndArrayLagrange = -(p.L1.*cos(theta1ArrayLagrange) + p.L2.*cos(theta2ArrayLagrange));

xOLagrange = zeros(length(tArrayLagrange), 1);
yOLagrange = zeros(length(tArrayLagrange), 1);

xEArrayLagrange = p.L1.*sin(theta1ArrayLagrange);
yEArrayLagrange = -(p.L1.*cos(theta1ArrayLagrange));

%The thing to actually plot
xPlotLagrange = [xOLagrange, xEArrayLagrange, xEndArrayLagrange]; 
yPlotLagrange = [yOLagrange, yEArrayLagrange, yEndArrayLagrange];

figure(5)
for i = 1:length(tArrayLagrange)
    plot(xPlotLagrange(i,:),yPlotLagrange(i,:),'LineWidth', 5)  
    axis equal
    axis ([-3 3 -3 3]);
    title('Animation of Double Pendulum-Lagrange')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end
figure(6);      %just the path of it
plot(xEndArrayLagrange, yEndArrayLagrange);
title('Position of the end of 2nd link over time-Lagrange')
xlabel('x Position')
ylabel('y Position')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

figure(7);
plot(tArrayLagrange, xEndArrayAMB - xEndArrayDAE);



function zDot = rhsAMB(z,p)
    %unpacks all the parameters into rhs function using the AMB method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta1dot = z(3); theta2dot = z(4);
    
    %getting the EOM from AMB, ie. the functions generated from 
    %dPendDeriveEOMamb file
    A = rhsStuffmassMatrixAMB(I1,I2,L1,d1,d2,m1,m2,theta1,theta2);
    b = rhsStuffbVectorAMB(L1,d1,d2,g,m1,m2,theta1,theta2,theta1dot,theta2dot);
    
    q = A\b;
    theta1doubledot = q(1);
    theta2doubledot = q(2);
    
    zDot = [theta1dot, theta2dot, theta1doubledot, theta2doubledot]'; 
end


function zDot = rhsDAE(z,p)
    %unpacks all the parameters into rhs function using the DAE method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta1dot = z(3); theta2dot = z(4);
    
    %getting the EOM from AMB, ie. dPendDeriveEOMdae file
    A = rhsStuffmassMatrixDAE(I1,I2,L1,d1,d2,m1,m2,theta1,theta2);
    b = rhsStuffbVectorDAE(L1,d1,d2,g,m1,m2,theta1,theta2,theta1dot,theta2dot);

    q = A\b;
    
    theta1doubledot =q(3);
    theta2doubledot =q(6);
    zDot = [theta1dot, theta2dot, theta1doubledot, theta2doubledot]'; 
end


function zDot = rhsLagrange(z,p)
    %unpacks all the parameters into rhs function using the Lagrange method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta1dot = z(3); theta2dot = z(4);
    
    %getting the EOM from AMB, ie. dPendDeriveEOMLagrange file
    M_L = rhsStuffmassMatrixLagrange(I1,I2,L1,d1,d2,m1,m2,theta1,theta2);
    b_L = rhsstuffbVectorLagrange(L1,d1,d2,g,m1,m2,theta1,theta2,theta1dot,theta2dot);
    
    q = M_L\b_L;
    theta1doubledot =q(1);
    theta2doubledot =q(2);
    
    zDot = [theta1dot, theta2dot, theta1doubledot, theta2doubledot]'; 
end

