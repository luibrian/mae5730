%Brian Lui
%bl569
%MAE5730
%HW#5

%Pendulum Numerics, Number 21
%Set up the pendulum in cartesian coordinates. Express the constant length 
%constraint as a set of linear equations restricting the acceleration. 
%Solve these (3 2nd order) DAE equations with numerical integration and 
%initial conditions and parameters of your choosing. No polar coordinates 
%allowed. Quantitatively compare your solution with a solution of the simple 
%pendulum equations (For the comparison you need to either compute x from 
%or vice versa. Integrate for a long enough time so you can detect drift 
%away from satisfying the kinematic constraint.

clear all;
close all;

p.m = 1;                %mass of the pendulum
p.g = 9.8;              %gravity

dur = 10;
npoints = 1001;
tspan = linspace(0,dur,npoints);

x0 = 1;
y0 = 1;
z0 = [x0;y0;0;0];         %[x,y,xDot, yDot]

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);

%This is weird because I have ex pointing down in my drawing, ey pointing
%to the right
xMatArray = zArray(:,2);
yMatArray = -zArray(:,1);

figure(1)
plot(xMatArray, yMatArray);
axis equal
title('Spatial path of pendulum- DAE')
xlabel('X position')
ylabel('Y Position')
hold on

% figure(4)
% for i = 1:length(tArray)
%     plot(xMatArray(i),yMatArray(i), '*');
%     hold on
%     axis([-2 2 -2 2]);
%     pause(0.001)
% end

%Solution of the simple pendulum equations
p.el = sqrt(x0^2 + y0^2);
theta0 = atan2(y0,x0);
z0ode = [theta0; 0];


g = @(t,z) rhsODE(z,p);
[tArrayODE, zArrayODE] = ode45(g, tspan, z0ode);

thetaArray = zArrayODE(:,1);
xODEArray = p.el .* sin(thetaArray);
yODEArray = -p.el .* cos(thetaArray);

figure(2);
plot(xODEArray, yODEArray);
title('Spatial path of pendulum-ode')
xlabel('X position')
ylabel('Y Position')
axis equal


xError = xMatArray - xODEArray;

figure(3)
plot(tArrayODE, xError);
title('Error in y position vs time')    %xposition, but y in plot
xlabel('Time')
ylabel('y Error')


function zDot = rhs(z,p)
    %DAE Method
    m = p.m; g = p.g; 
    x = z(1); y = z(2); xDot = z(3); yDot = z(4);
    
    el0 = sqrt(x^2 + y^2);
    
    
    %Do Av = b, where v = [xDoubleDot;yDoubleDot;T]
    %T = tension
    A = [m, 0, x/el0;
        0, m, y/el0; 
        x, y, 0];
    b = [m*g; 0; -(xDot^2+yDot^2)];
    
    v = A\b;
    xDoubleDot = v(1);
    yDoubleDot = v(2);
    
    zDot = [xDot; yDot; xDoubleDot; yDoubleDot];

end


function zDot = rhsODE(z,p)
    m = p.m; g=p.g; el = p.el;
    theta = z(1); thetaDot = z(2);
    
    thetaDoubleDot = -g/el *sin(theta);
    
    zDot = [thetaDot; thetaDoubleDot];
end
