%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%4 Bar Linkage. Hinges at origin O, then elbow E, then elbow F. Constrained
%at the other end at point H. g=1.
%Neglect all friction and assume there are no joint motors. 

%Parameters
p.L1 = 1; p.L2 = 1; p.L3 = 2;
p.d1 = 0.5; p.d2 = 0.5; p.d3 = 0.5;
p.m1 = 1; p.m2 = 1; p.m3 = 1;
%p.I1 = 1/12*p.m1*p.d1^2; p.I2 = 1/12*p.m2*p.d2^2; 
p.g = 1; p.I1 = 0.2; p.I2 = 0.3; p.I3 = 0.4;

%setting up tspan
dur = 100;
npoints = 200;
tspan = linspace(0, dur, npoints);

%initial condition
%They must satisfy constraint
%z0 = [pi/4; pi; pi;0;0;0];      %stationary condition
theta1_0 = pi - pi/30;
theta2_0 = pi/2;
theta3_0 = 0;
theta1dot_0 = 0;
theta2dot_0 = 0;
theta3dot_0 = 0;

%{
x1_0 = p.d1*cos(theta1_0);
y1_0 = p.d1*sin(theta1_0);
x2_0 = p.L1*cos(theta1_0) + p.d2*cos(theta2_0);
y2_0 = p.L1*sin(theta1_0) + p.d2*sin(theta2_0);
x3_0 = p.L1*cos(theta1_0) + p.L2*cos(theta2_0) + p.d3*cos(theta3_0);
y3_0 = p.L1*sin(theta1_0) + p.L2*cos(theta2_0) + p.d3*cos(theta3_0);

x1dot_0 = -p.d1*theta1dot_0*sin(theta1_0);
y1dot_0 = p.d1*theta1dot_0*cos(theta1_0);
x2dot_0 = -p.L1*theta1dot_0*sin(theta1_0) -p.d2*theta2dot_0*sin(theta2_0);
y2dot_0 = p.L1*theta1dot_0*cos(theta1_0) +p.d2*theta2dot_0*cos(theta2_0);
x3dot_0 = -p.L1*theta1dot_0*sin(theta1_0) -p.L2*theta2dot_0*sin(theta2_0) -p.d3*theta3dot_0*sin(theta3_0); 
y3dot_0 = p.L1*theta1dot_0*cos(theta1_0) +p.L2*theta2dot_0*cos(theta2_0) + p.d3*theta3dot_0*cos(theta3_0);

pos = [x1_0; y1_0; x2_0; y2_0; x3_0; y3_0];
vel = [x1dot_0; y1dot_0; x2dot_0; y2dot_0; x3dot_0; y3dot_0];

z0 = [theta1_0; theta2_0; theta3_0; theta1dot_0; theta2dot_0; theta3dot_0;...
    pos; vel];      
%}
z0 = [theta1_0; theta2_0; theta3_0; theta1dot_0; theta2dot_0; theta3dot_0];

%Setting tolerance
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);


%Solving and plotting
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
 
%
figure(1)
for i = 1:length(tArrayDAE)
    plot(xPlotDAE(i,:),yPlotDAE(i,:),'LineWidth', 5)  
    axis equal
    axis ([-4 4 -4 4]);
    title('Animation of Triple Pendulum-DAE')
    xlabel('x Position')
    ylabel('y Position')
    pause(0.00001);
    shg;
end

%{
%Energy check
x1_Array = zArrayDAE(:,7);
y1_Array = zArrayDAE(:,8);
x2_Array = zArrayDAE(:,9);
y2_Array = zArrayDAE(:,10);
x3_Array = zArrayDAE(:,11);
y3_Array = zArrayDAE(:,12);
x1dot_Array = zArrayDAE(:,13);
y1dot_Array = zArrayDAE(:,14);
x2dot_Array = zArrayDAE(:,15);
y2dot_Array = zArrayDAE(:,16);
x3dot_Array = zArrayDAE(:,17);
y3dot_Array = zArrayDAE(:,18);

%pot1 = p.m1.*p.g.*(p.d1.*cos(theta1ArrayDAE));
%pot2 = p.m2.*p.g.*(p.d2.*cos(theta2ArrayDAE) .+ p.L1.*cos(theta1ArrayDAE));
%pot3 = p.m3.*p.g.*(p.d3.*cos(theta3ArrayDAE) .+ p.L2.*cos(theta2ArrayDAE) + p.L1*cos(theta1ArrayDAE));
pot1 = p.m1.*p.g.*x1_Array;
pot2 = p.m2.*p.g.*x2_Array;
pot3 = p.m3.*p.g.*x3_Array;
potTot = pot1+pot2+pot3;
vel1 = 1/2.* p.m1.* (x1dot_Array + y1dot_Array).^2 + 1/2.*p.I1.*theta1ArrayDAE.^2;
vel2 = 1/2.* p.m2.* (x2dot_Array + y2dot_Array).^2 + 1/2.*p.I2.*theta2ArrayDAE.^2;
vel3 = 1/2.* p.m3.* (x3dot_Array + y3dot_Array).^2 + 1/2.*p.I3.*theta3ArrayDAE.^2;
velTot = vel1 + vel2 + vel3; 
Etot = potTot + velTot;

figure(5);
plot(tArrayDAE, Etot);
%}


function zDot = rhsDAE(z,p)
    %unpacks all the parameters into rhs function using the DAE method
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end
    
    %unpack the state
    theta1 = z(1); theta2 = z(2); theta3 = z(3);
    theta1dot = z(4); theta2dot = z(5); theta3dot = z(6);
    
%    velXY = z(13:18)';
    %getting the EOM from AMB, ie. dPendDeriveEOMdae file
    A = fourBar_rhsStuffmassMatrixDAE(I1,I2,I3,L1,L2,L3,d1,d2,d3,m1,m2,m3,theta1,theta2,theta3);
    b = fourBar_rhsStuffbVectorDAE(L1,L2,L3,d1,d2,d3,g,m1,m2,m3,theta1,theta2,theta3,theta1dot,theta2dot,theta3dot);
    
    q = A\b;
    
    
    theta1doubledot =q(3);
    theta2doubledot =q(6);
    theta3doubledot =q(9);

    
%    accelXY = [q(1) q(2) q(4) q(5) q(7) q(8)];
    zDot = [theta1dot, theta2dot, theta3dot,...
        theta1doubledot, theta2doubledot theta3doubledot]';
%    zDot = [theta1dot, theta2dot, theta3dot,...
%        theta1doubledot, theta2doubledot theta3doubledot,...
%        velXY, accelXY]'; 
end

