%Brian Lui
%Bl569
%MAE5780 Intermediate Dynamics
%9/19/18
%HW3

clear all
close all

%#14- Mechanics of 2(+) particles

%Part a/b- 2 unequal particles with mass m1 and m2, find period of circular
%motion if the distance btwn them is d and the force between them is F =
%Gm1m2/r^2.

%adjustable paramters
p.G = 1;               %gravitational constant
p.m1 = 1;               %mass 1
p.m2 = 2;               %mass 2
p.d = 5;                %distance btwn the m1 and m2
p.y = p.m2/(p.m1+p.m2)*p.d;     %dist from m1 to CM


theoreticalPeriod = 2*pi*sqrt(p.d^3/(p.G*p.m1*p.m2));

%must choose a value in which the circles are complete
dur = theoreticalPeriod;  
npoints = 1001;
tspanA = linspace(0,dur,npoints);

%must choose initial starting points that are pi rads apart for putting CM
%about the the origin
%This puts the CM at the origin. The distance of the second circle from the
%CM of the system = d-r1 (distance btwn the two points - distance of 1
%circle from the CM)
r1 = p.y;
r2 = p.d-r1;
z0A = [r1; 0; r2; pi];      %r1,theta1, r2, theta2


%Doing ode45 and extracting its relevant elements
f = @(t,z) rhsA(z,p);
[tArrayA, zArrayA] = ode45(f,tspanA, z0A);
r1ArrayA = zArrayA(:,1);
theta1ArrayA = zArrayA(:,2);
r2ArrayA = zArrayA(:,3);
theta2ArrayA = zArrayA(:,4);

x1A = zeros(length(tArrayA),1);
y1A = zeros(length(tArrayA),1);
x2A = zeros(length(tArrayA),1);
y2A = zeros(length(tArrayA),1);


for i = 1:length(tArrayA)
    x1A(i) = r1ArrayA(i)*cos(theta1ArrayA(i));
    y1A(i) = r1ArrayA(i)*sin(theta1ArrayA(i));
    x2A(i) = r2ArrayA(i)*cos(theta2ArrayA(i));
    y2A(i) = r2ArrayA(i)*sin(theta2ArrayA(i));
end

figure(1);
plot(x1A(1:end),y1A(1:end));
hold on;
plot(x2A,y2A);
axis([-5 5 -5 5]);

%seeing how it gets plotted
figure(2);
for i=1:length(tArrayA)
    plot(x1A(i), y1A(i), 'g*');
    hold on;
    plot(x2A(i), y2A(i), 'b*');
    axis([-5 5 -5 5]);
    pause(0.000001);
end



function zDot = rhsA(z, p)
    %Done in polar
    r1 = z(1); theta1 = z(2); r2 = z(3); theta2 =z(4);
    
    G = p.G; m1 = p.m1; m2 = p.m2; d = p.d; y =p.y; %parameters
    
    
    r1Dot = 0;  %circular motion, cannot have change in radius
    r2Dot = 0;
    
    theta1Dot = sqrt(G*(m1+m2)/d^3);
    theta2Dot = sqrt(G*m1/(d^2 *(d- y)));
    
    zDot = [r1Dot; theta1Dot; r2Dot; theta2Dot];
end