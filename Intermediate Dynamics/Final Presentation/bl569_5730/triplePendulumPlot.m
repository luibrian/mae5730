%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%Triple pendulum. Hinges at origin O, then elbow E, then elbow F.
%Neglect all friction and assume there are no joint motors. 

%Parameters
p.L1 = 1; p.L2 = 1; p.L3 = 1;
p.d1 = 0.5; p.d2 = 0.5; p.d3 = 0.5;
p.m1 = 1; p.m2 = 1; p.m3 = 1;
%p.I1 = 1/12*p.m1*p.d1^2; p.I2 = 1/12*p.m2*p.d2^2; 
p.g = 1; p.I1 = 0.2; p.I2 = 0.3; p.I3 = 0.4;

%setting up tspan
dur = 20;
npoints = 200;
tspan = linspace(0, dur, npoints);

%initial condition
%z0 = [pi/2; pi; pi; 0; 0; 0];      
%z0 = [pi/4; pi/2; 3*pi/4;0;0;0];
%z0 = [pi/2; pi; pi;0;0;0];
z0 = [pi; pi/3; pi/3;0;0;0];

%Setting tolerance
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the NE Method
fAMB = @(t,z) rhsAMB(z,p);
[tArrayAMB, zArrayAMB] = ode45(fAMB, tspan, z0, options);

theta1ArrayAMB = zArrayAMB(:,1);
theta2ArrayAMB = zArrayAMB(:,2);
theta3ArrayAMB = zArrayAMB(:,3);

%Getting the points of interest
xEArrayAMB = p.L1.*sin(theta1ArrayAMB);
yEArrayAMB = -(p.L1.*cos(theta1ArrayAMB));

xFArrayAMB = p.L1.*sin(theta1ArrayAMB) + p.L2.*sin(theta2ArrayAMB);
yFArrayAMB = -(p.L1.*cos(theta1ArrayAMB) + p.L2.*cos(theta2ArrayAMB));

xEndArrayAMB = xFArrayAMB + p.L3.*sin(theta3ArrayAMB);
yEndArrayAMB = yFArrayAMB - p.L3.*cos(theta3ArrayAMB);

xOAMB = zeros(length(tArrayAMB), 1);
yOAMB = zeros(length(tArrayAMB), 1);

%The thing to actually plot
xPlotAMB = [xOAMB, xEArrayAMB, xFArrayAMB, xEndArrayAMB]; 
yPlotAMB = [yOAMB, yEArrayAMB, yFArrayAMB, yEndArrayAMB];
 
%
figure(1)
for i = 1:length(tArrayAMB)
    plot(xPlotAMB(i,:),yPlotAMB(i,:),'LineWidth', 5)  
    axis equal
    axis ([-6 6 -6 6]);
    title('Animation of Triple Pendulum-AMB')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end
%
%
figure(2);      %just the path of it
plot(xEndArrayAMB, yEndArrayAMB);
%title('Position of the end of 3rd link over time-AMB')
xlabel('x Position')
ylabel('y Position')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the DAE Method
fDAE = @(t,z) rhsDAE(z,p);
[tArrayDAE, zArrayDAE] = ode45(fDAE, tspan, z0, options);

theta1ArrayDAE = zArrayDAE(:,1);
theta2ArrayDAE = zArrayDAE(:,2);
theta3ArrayDAE = zArrayDAE(:,3);

%Getting the points of interest
xEArrayDAE = p.L1.*sin(theta1ArrayDAE);
yEArrayDAE = -(p.L1.*cos(theta1ArrayDAE));

%xFArrayDAE = p.L1.*sin(theta1ArrayDAE) + p.L2.*sin(theta2ArrayDAE);
%yFArrayDAE = -(p.L1.*cos(theta1ArrayDAE) + p.L2.*cos(theta2ArrayDAE));
xFArrayDAE = xEArrayDAE + p.L2.*sin(theta2ArrayDAE);
yFArrayDAE = yEArrayDAE - p.L2.*cos(theta2ArrayDAE);

xEndArrayDAE = xFArrayDAE + p.L3.*sin(theta3ArrayDAE);
yEndArrayDAE = yFArrayDAE - p.L3.*cos(theta3ArrayDAE);

xODAE = zeros(length(tArrayDAE), 1);
yODAE = zeros(length(tArrayDAE), 1);


%The thing to actually plot
xPlotDAE = [xODAE, xEArrayDAE, xFArrayDAE, xEndArrayDAE]; 
yPlotDAE = [yODAE, yEArrayDAE, yFArrayDAE, yEndArrayDAE];
 
%{
figure(3)
for i = 1:length(tArrayDAE)
    plot(xPlotDAE(i,:),yPlotDAE(i,:),'LineWidth', 5)  
    axis equal
    axis ([-6 6 -6 6]);
    title('Animation of Triple Pendulum-DAE')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end
%}
figure(4);      %just the path of it
plot(xEndArrayDAE, yEndArrayDAE);
title('Position of the end of 3rd link over time-DAE')
xlabel('x Position')
ylabel('y Position')
%}

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving using the Lagrange Method
fLagrange = @(t,z) rhsLagrange(z,p);
[tArrayLagrange, zArrayLagrange] = ode45(fLagrange, tspan, z0, options);

theta1ArrayLagrange = zArrayLagrange(:,1);
theta2ArrayLagrange = zArrayLagrange(:,2);
theta3ArrayLagrange = zArrayLagrange(:,3);

%Getting the points of interest
xFArrayLagrange = p.L1.*sin(theta1ArrayLagrange) + p.L2.*sin(theta2ArrayLagrange);
yFArrayLagrange = -(p.L1.*cos(theta1ArrayLagrange) + p.L2.*cos(theta2ArrayLagrange));

xEArrayLagrange = p.L1.*sin(theta1ArrayLagrange);
yEArrayLagrange = -(p.L1.*cos(theta1ArrayLagrange));

xEndArrayLagrange = xFArrayLagrange + p.L3.*sin(theta3ArrayLagrange);
yEndArrayLagrange = yFArrayLagrange - p.L3.*cos(theta3ArrayLagrange);


xOLagrange = zeros(length(tArrayLagrange), 1);
yOLagrange = zeros(length(tArrayLagrange), 1);


%The thing to actually plot
xPlotLagrange = [xOLagrange, xEArrayLagrange, xFArrayLagrange, xEndArrayLagrange]; 
yPlotLagrange = [yOLagrange, yEArrayLagrange, yFArrayLagrange, yEndArrayLagrange];
%{
figure(5)
for i = 1:length(tArrayLagrange)
    plot(xPlotLagrange(i,:),yPlotLagrange(i,:),'LineWidth', 5)  
    axis equal
    axis ([-6 6 -6 6]);
    title('Animation of Triple Pendulum-Lagrange')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;    
end
%}
figure(6);      %just the path of it
plot(xEndArrayLagrange, yEndArrayLagrange);
title('Position of the end of 3rd link over time-Lagrange')
xlabel('x Position')
ylabel('y Position')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%
figure(7);
plot(tArrayLagrange, xEndArrayAMB- xEndArrayLagrange)
xlabel('Time(s)');
ylabel('x_{NE} - x_{Lagrange}');
%
%
figure(8);
plot(tArrayAMB, xEndArrayAMB - xEndArrayDAE);
xlabel('Time(s)');
ylabel('x_{NE} - x_{DAE}');
%


function zDot = rhsAMB(z,p)
    %unpacks all the parameters into rhs function using the NE method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta3 = z(3);
    theta1dot = z(4); theta2dot = z(5); theta3dot = z(6);
    %W = z(7);   %For energy check
    
    
    %getting the EOM from AMB, ie. the functions generated from 
    %dPendDeriveEOMamb file
    A = rhsStuffmassMatrixAMB(I1,I2,I3,L1,L2,d1,d2,d3,m1,m2,m3,theta1,theta2,theta3);
    b = rhsStuffbVectorAMB(L1,L2,d1,d2,d3,g,m1,m2,m3,theta1,theta2,theta3,theta1dot,theta2dot,theta3dot);    
    q = A\b;
    theta1doubledot = q(1);
    theta2doubledot = q(2);
    theta3doubledot = q(3);
    
    %Energy Check
    %Wdot = P;
    
    zDot = [theta1dot, theta2dot, theta3dot,...
        theta1doubledot, theta2doubledot theta3doubledot]'; 
end

function zDot = rhsDAE(z,p)
    %unpacks all the parameters into rhs function using the DAE method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta3 = z(3);
    theta1dot = z(4); theta2dot = z(5); theta3dot = z(6);
    
    %getting the EOM from AMB, ie. dPendDeriveEOMdae file
    A = rhsStuffmassMatrixDAE(I1,I2,I3,L1,L2,d1,d2,d3,m1,m2,m3,theta1,theta2,theta3);
    b = rhsStuffbVectorDAE(L1,L2,d1,d2,d3,g,m1,m2,m3,theta1,theta2,theta3,theta1dot,theta2dot,theta3dot);

    q = A\b;
    
    theta1doubledot =q(3);
    theta2doubledot =q(6);
    theta3doubledot =q(9);
    zDot = [theta1dot, theta2dot, theta3dot,...
        theta1doubledot, theta2doubledot theta3doubledot]'; 
end


function zDot = rhsLagrange(z,p)
    %unpacks all the parameters into rhs function using the Lagrange method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta3 = z(3);
    theta1dot = z(4); theta2dot = z(5); theta3dot = z(6);
    
    
    %getting the EOM from AMB, ie. dPendDeriveEOMLagrange file
    M_L = rhsStuffmassMatrixLagrange(I1,I2,I3,L1,L2,d1,d2,d3,m1,m2,m3,theta1,theta2,theta3);
    b_L = rhsstuffbVectorLagrange(L1,L2,d1,d2,d3,g,m1,m2,m3,theta1,theta2,theta3,theta1dot,theta2dot,theta3dot);
    
    q = M_L\b_L;
    theta1doubledot =q(1);
    theta2doubledot =q(2);
    theta3doubledot =q(3);
    
    zDot = [theta1dot, theta2dot, theta3dot,...
        theta1doubledot, theta2doubledot theta3doubledot]'; 
end

