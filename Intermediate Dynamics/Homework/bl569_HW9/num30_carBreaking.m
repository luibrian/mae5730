%Brian Lui
%MAE 5730
%HW9
%10/27/18

clear all
close all

%Breaking Stability 2D. 
%Assume steering is locked (either front brakes or rear breaks) and
%straight ahead. For simplicity, assume center of mass is at ground heigh
%between the front and back wheels. Assume that the locked wheels act as a
%single dragging point on the centerline of the car midway between the
%wheels. 

%Rear wheels breaking
%Parameters
p.m = 1;                    %mass of car
p.I = 1;                    %moment of inertia of car
p.d1 = 0.5;                 %distance between rear wheels and CM
p.d2 = 0.3;                 %distance between front wheels and CM
p.muu = 1;                  %Coefficient of friction
p.L = p.d1 + p.d2;          %total distance between the 
p.g = 1;

%Timespan
dur = 5;
npoints = 301;
tspan = linspace(0,dur,npoints);

%Initial Conditions
z0 = [0;0;0;5;1]; %x;y;theta;v;omega
%For rear brakes applied, x and y are the location of the front wheels
%For front brakes applied, x and y are the location of the rear wheels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rear brakes
fRear = @(t,z) rhsRear(z,p);
[tArrayRear, zArrayRear] = ode45(fRear,tspan, z0);

%Extracting the relevant vectors
xArrayRearA = zArrayRear(:,1);
yArrayRearA = zArrayRear(:,2);
thetaArrayRear = zArrayRear(:,3);
vArrayRear = zArrayRear(:,4);
omegaArrayRear = zArrayRear(:,5);


xArrayRearG = xArrayRearA - p.d2.* cos(thetaArrayRear);
xArrayRearC = xArrayRearA - p.L.* cos(thetaArrayRear);
yArrayRearG = yArrayRearA - p.d2.* sin(thetaArrayRear);
yArrayRearC = yArrayRearA - p.L.* sin(thetaArrayRear);
    

%Uncomment me to see the path 
%
figure(4);
%Animating path of car
for i=1:length(tArrayRear)
    
    plot(xArrayRearA(i), yArrayRearA(i),'r*');
    hold on;
    plot(xArrayRearG(i), yArrayRearG(i),'g*');
    plot(xArrayRearC(i), yArrayRearC(i),'b*');

    axis equal;
    axis ([-1 4 -1 4]);
    title('Rear brakes, Position in space as t increases');
    pause(0.01);
    hold off;
end
%

%plotting things for rear brakes
figure(1);
%Plotting path of car
plot(xArrayRearA, yArrayRearA,'r');
hold on;
plot(xArrayRearG, yArrayRearG,'g');
plot(xArrayRearC, yArrayRearC,'b');

axis equal;
%axis ([-1 4 -1 4]);

%Starting points of car
plot(xArrayRearA(1), yArrayRearA(1), 'r*');
plot(xArrayRearG(1), yArrayRearG(1), 'g*');
plot(xArrayRearC(1), yArrayRearC(1), 'b*');
%Ending points of car
plot(xArrayRearA(end), yArrayRearA(end), 'r*');
plot(xArrayRearG(end), yArrayRearG(end), 'g*');
plot(xArrayRearC(end), yArrayRearC(end), 'b*');
axis equal;

%Labelling plot
title('Rear brakes, Position in space as t increases');
xlabel('X Position');
ylabel('Y Position');
legend('Front wheels','CM of car','Rear wheels', 'Location', 'NW');

%{
figure(2);
subplot(2,1,1);
plot(tArrayRear, vArrayRear);
title('Rear brakes, velocity vs time');
xlabel('Time');
ylabel('Velocity of skating wheel');
subplot(2,1,2);
plot(tArrayRear, omegaArrayRear);
title('Rear brakes, angular velocity vs time');
xlabel('Time');
ylabel('\omega of car');

%Plot of the distance between front and back wheels vs time 
figure(3);
plot(tArrayRear, sqrt((xArrayRearA- xArrayRearC).^2 + (yArrayRearA- yArrayRearC).^2));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Front brakes
fFront = @(t,z) rhsFront(z,p);
[tArrayFront, zArrayFront] = ode45(fFront,tspan, z0);

xArrayFrontA    = zArrayFront(:,1);
yArrayFrontA    = zArrayFront(:,2);
thetaArrayFront = zArrayFront(:,3);

xArrayFrontG = xArrayFrontA + p.d1.* cos(thetaArrayFront);
xArrayFrontC = xArrayFrontA + p.L.* cos(thetaArrayFront);
yArrayFrontG = yArrayFrontA + p.d1.* sin(thetaArrayFront);
yArrayFrontC = yArrayFrontA + p.L.* sin(thetaArrayFront);

vArrayFront = zArrayFront(:,4);
omegaArrayFront = zArrayFront(:,5);

%plotting things for front brakes
figure(5);
%Plotting path of car
plot(xArrayFrontA, yArrayFrontA,'r');
hold on;
plot(xArrayFrontG, yArrayFrontG,'g');
plot(xArrayFrontC, yArrayFrontC,'b');
%Starting points of car
plot(xArrayFrontA(1), yArrayFrontA(1), 'r*');
plot(xArrayFrontG(1), yArrayFrontG(1), 'g*');
plot(xArrayFrontC(1), yArrayFrontC(1), 'b*');
%Ending points of car
plot(xArrayFrontA(end), yArrayFrontA(end), 'r*');
plot(xArrayFrontG(end), yArrayFrontG(end), 'g*');
plot(xArrayFrontC(end), yArrayFrontC(end), 'b*');
axis equal;
axis ([-1 15 -1 15]);

%Labelling plot
title('Front brakes, Position in space as t increases');
xlabel('X Position');
ylabel('Y Position');
legend('Rear wheels','CM of car','Front wheels', 'Location', 'NW');
%{
figure(6);
subplot(2,1,1);
plot(tArrayFront, vArrayFront);
title('Front wheels, velocity vs time');
xlabel('Time');
ylabel('Velocity of skating wheel');
subplot(2,1,2);
plot(tArrayFront, omegaArrayFront);
title('Front wheels, angular velocity vs time');
xlabel('Time');
ylabel('\omega of car');
%}
%Uncomment me to see the path 
%{
figure(7);
%Animating path of car
for i=1:length(tArrayFront)
    
    plot(xArrayFrontA(i), yArrayFrontA(i),'r*');
    hold on;
    plot(xArrayFrontG(i), yArrayFrontG(i),'g*');
    plot(xArrayFrontC(i), yArrayFrontC(i),'b*');

    axis equal;
    axis ([-1 15 -1 15]);
    
    pause(0.01);
    hold off;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function zdot = rhsRear(z,p)
    %unpacks all the parameters into rhs function using the DAE method
    %Rear brakes locked, front wheels skate. Ie. Rear wheel brakes.
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end    
    
    x = z(1); 
    y = z(2);
    theta = z(3);
    v = z(4);           %velocity of car in lambda direction
    omega = z(5);       %angular velocity of car
    
    %by definition
    xDot = v*cos(theta);
    yDot = v*sin(theta);
    thetaDot = omega;
    
    %EOMs
    totalVmag = sqrt(v^2 + L^2 *omega^2);
    vDot = -(muu*g*v/totalVmag + d2*omega^2);
    omegaDot = 1/(m*d2^2 + I) * (m*d2*(v*omega) - muu*m*g*L^2*omega/totalVmag);
    
    zdot = [xDot; yDot; thetaDot; vDot; omegaDot];
end


function zdot = rhsFront(z,p)
    %unpacks all the parameters into rhs function using the DAE method
    %Rear brakes skate, front wheels lock. Ie. Front wheel brakes
    names = fieldnames(p);
    for i = 1:length(names)
        eval([names{i} '= p.' names{i} ';']);
    end    
    
    x = z(1); 
    y = z(2);
    theta = z(3);
    v = z(4);           %velocity of car in lambda direction
    omega = z(5);       %angular velocity of car
    
    %by definition
    xDot = v*cos(theta);
    yDot = v*sin(theta);
    thetaDot = omega;
    
    %EOMs
    totalVmag = sqrt(v^2 + L^2 *omega^2);
    vDot = -(muu*g*v/totalVmag) + d1*omega^2;
    omegaDot = -1/(m*d1^2 + I)*((m*d1*v*omega) +...
        muu*m*g*L^2*omega/totalVmag);
    
    zdot = [xDot; yDot; thetaDot; vDot; omegaDot];
end
