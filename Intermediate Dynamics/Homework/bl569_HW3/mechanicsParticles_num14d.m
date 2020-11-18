%Brian Lui
%Bl569
%MAE5780 Intermediate Dynamics
%9/19/18
%HW3

clear all
close all

%#14d- Mechanics of 2(+) particles
%Part c/d- 3 equal particles m1=m2=m3=1, G=10, what is the angular speed
%for circular motion on a circle with diameter of d=3?

%adjustable paramters
p.G = 10;               %gravitational constant
p.m1 = 1;               %mass 1
p.m2 = 1;               %mass 2
p.m3 = 1;               %mass 3
p.d = 3;                %diameter of the circle

%theoreticalPeriod = 2*pi*sqrt(p.d^3/(p.G*p.m1*p.m2));
%Constant rotation throughout
thetaDot = sqrt(8*p.G*p.m1/(sqrt(3)*p.d^3));

period = 2*pi/thetaDot;

%must choose a value in which the circles are complete
%This should go exactly so that 1 circle is filled in
dur = period/3;  
npoints = 1001;
tspan = linspace(0,dur,npoints);

%must choose initial starting points that are on equidistant points on
%circle
%This puts the CM at the origin. 

r = p.d/2;              %radius of the circle

%initial velocity conditions at mass 1
v1_0 = r*thetaDot;
v1_0x = -v1_0;
v1_0y = 0;

%Getting the initial velocity at mass 2 and 3   
rotationMatrixTo2 = [cosd(120), -sind(120); sind(120), cosd(120)];
rotationMatrixTo3 = [cosd(-120), -sind(-120); sind(-120), cosd(-120)];
v2 = rotationMatrixTo2 * [v1_0x; v1_0y];
v3 = rotationMatrixTo3 * [v1_0x; v1_0y];
v2_0x = v2(1);
v2_0y = v2(2);
v3_0x = v3(1);
v3_0y = v3(2);



%[x1;y1;x2;y2;x3;y3;x1Dot;y1Dot;x2Dot;y2Dot;x3Dot;y3Dot]
p0 = [0; r;...                      %x1,y1
    -r*cosd(30); -r*sind(30);...    %x2,y2
    r*cosd(30) ; -r*sind(30)];      %x3,y3
p0Dot = [v1_0x; v1_0y;              %xDot1, yDot1
    v2_0x; v2_0y;                   %xDot2, yDot2
    v3_0x; v3_0y];                  %x3Dot, y3Dot
z0 = [p0; p0Dot];
    

%Doing ode45 and extracting its relevant elements
options = odeset('Abstol',1e-8,'Reltol',1e-8);

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f,tspan,z0,options);
x1Array = zArray(:,1);
y1Array = zArray(:,2);
x2Array = zArray(:,3);
y2Array = zArray(:,4);
x3Array = zArray(:,5);
y3Array = zArray(:,6);

%Plotting
figure(1);
plot(x1Array,y1Array,'r');
hold on;
plot(x2Array,y2Array,'b');
plot(x3Array,y3Array,'g');

% for i=1:length(tArray)
%     plot(x1Array(i),y1Array(i),'g*');
%     hold on;
%     plot(x2Array(i),y2Array(i),'b*');
%     plot(x3Array(i),y3Array(i),'r*');
%     pause(0.00001)
% end

%Error
x1err = x1Array(end)-p0(3);
y1err = y1Array(end)-p0(4);
x2err = x2Array(end)-p0(5);
y2err = y2Array(end)-p0(6);
x3err = x3Array(end)-p0(1);
y3err = y3Array(end)-p0(2);

[x1err; y1err; x2err; y2err; x3err; y3err]



function zDot = rhs(z, p)
    %Done in cartesian- states
    x1 = z(1); y1 = z(2); 
    x2 = z(3); y2 =z(4); 
    x3 = z(5); y3 = z(6);
    x1Dot = z(7); y1Dot = z(8); 
    x2Dot = z(9); y2Dot = z(10); 
    x3Dot = z(11); y3Dot = z(12);
    
    %reading in parameters
    G = p.G; 
    m1 = p.m1; m2 = p.m2; m3 = p.m3;
    
    %distances (these are honestly all equal)
    el_12 = sqrt((y2-y1)^2 + (x2-x1)^2);
    el_13 = sqrt((y3-y1)^2 + (x3-x1)^2);
    el_23 = sqrt((y3-y2)^2 + (x3-x2)^2);

    %equations of motions
    x1DoubleDot = G*m2/(el_12)^3 *(x2-x1) +  G*m3/(el_13)^3 *(x3-x1);
    y1DoubleDot = G*m2/(el_12)^3 *(y2-y1) +  G*m3/(el_13)^3 *(y3-y1);
    x2DoubleDot = G*m1/(el_12)^3 *(x1-x2) +  G*m3/(el_23)^3 *(x3-x2);
    y2DoubleDot = G*m1/(el_12)^3 *(y1-y2) +  G*m3/(el_23)^3 *(y3-y2);
    x3DoubleDot = G*m1/(el_13)^3 *(x1-x3) +  G*m2/(el_23)^3 *(x2-x3);
    y3DoubleDot = G*m1/(el_13)^3 *(y1-y3) +  G*m2/(el_23)^3 *(y2-y3);

    %returning the zDot
    zDot = [x1Dot; y1Dot; x2Dot; y2Dot; x3Dot; y3Dot;...
        x1DoubleDot; y1DoubleDot;...
        x2DoubleDot; y2DoubleDot;...
        x3DoubleDot; y3DoubleDot];
end