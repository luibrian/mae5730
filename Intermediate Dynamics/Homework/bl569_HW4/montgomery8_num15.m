%Brian Lui
%MAE 5730
%HW4
%9/26/18

clear all
close all

%Initial Conditions
x1 = -0.97000436;
y1 = 0.24308753;
x2 = -x1;
y2 = -y1;
x3 = 0;
y3 = 0;

vx3 = 0.93240737;
vy3 = 0.86473146;
vx1 = -vx3/2;
vy1 = -vy3/2;
vx2 = -vx3/2;
vy2 = -vy3/2;

z0 = [x1;y1;x2;y2;x3;y3; vx1;vy1;vx2;vy2;vx3;vy3];

%parameters
p.m1 = 1; p.m2 = 1; p.m3 = 1; p.G = 1;


durArray = [2.1, 10, 10, 100];
npoints = 1001;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
for i = 1:length(durArray)
    if i>2
        z0(1) = z0(1)-0.1;
        z0(2) = z0(2)-0.1;
    end
    dur = durArray(i);
    tspan = linspace(0,dur, npoints);
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan, z0, options);

    x1Array = zArray(:,1);
    y1Array = zArray(:,2);
    x2Array = zArray(:,3);
    y2Array = zArray(:,4);
    x3Array = zArray(:,5);
    y3Array = zArray(:,6);

    figure(i);
    plot(x1Array,y1Array, 'b');
    hold on;
    plot(x2Array,y2Array, 'g');
    plot(x3Array,y3Array, 'r--');
end


function zDot = rhs(z,p)
    m1 = p.m1; m2 = p.m2; m3 = p.m3; G = p.G;
    
    x1 = z(1); y1 = z(2); x2 = z(3); y2 = z(4); x3 = z(5); y3 = z(6);
    x1Dot = z(7); y1Dot = z(8); 
    x2Dot = z(9); y2Dot = z(10); 
    x3Dot = z(11); y3Dot = z(12);
    
    el_12 = sqrt((x2-x1)^2 + (y2-y1)^2);
    el_13 = sqrt((x3-x1)^2 + (y3-y1)^2);
    el_23 = sqrt((x3-x2)^2 + (y3-y2)^2);
    
    %equations of motions
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