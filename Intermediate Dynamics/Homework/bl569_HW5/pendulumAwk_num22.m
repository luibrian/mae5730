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

p.g = 9.8;              %gravity
p.el = 1;               %length of pendulum

dur = 10;
npoints = 1001;
tspan = linspace(0,dur,npoints);

y0 = 0.5;
y0Dot = 0;
z0 = [y0;y0Dot];         %[y; yDot]

f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);

%This is weird because I have ex pointing down in my drawing, ey pointing
%to the right
xMatArray = zArray(:,1);
yMatArray = -sqrt(p.el^2 - zArray(:,1).^2);
figure(1)
plot(xMatArray, yMatArray);
title('Spatial path of pendulum- DAE')
xlabel('X position')
ylabel('Y Position')
axis equal
hold on


%%%Animating the pendulum
% figure(4)
% for i = 1:length(tArray)
%     plot(xMatArray(i),yMatArray(i), '*');
%     hold on
%     axis([-2 2 -2 2]);
%     pause(0.001)
% end

%Solution of the simple pendulum equations
x0 = sqrt(p.el^2 - y0^2);
theta0 = atan2(y0,x0);
z0ode = [theta0; 0];

g = @(t,z) rhsODE(z,p);
[tArrayODE, zArrayODE] = ode45(g, tspan, z0ode);

thetaArray = zArrayODE(:,1);
xODEArray = p.el .* sin(thetaArray);
yODEArray = -p.el .* cos(thetaArray);

figure(2);
plot(xODEArray, yODEArray);
title('Spatial path of pendulum- ODE')
xlabel('X position')
ylabel('Y Position')
axis equal


%checking out the error in y position- comparing the two methods
yError = yMatArray - yODEArray;
figure(3)
plot(tArrayODE, yError);
title('Error in x position vs time')
xlabel('Time')
ylabel('x Error')


function zDot = rhs(z,p)
    %DAE Method- 2nd order equation using the horizontal position y as
    %state
    g=p.g; el = p.el;
    y = z(1); yDot = z(2);
    
    x = sqrt(el^2 - y^2);
    yDoubleDot = -y*(el^2/x^3 * yDot^2 + g) / (x+ y^2/x);
    
    zDot = [yDot; yDoubleDot];
end



function zDot = rhsODE(z,p)
    g=p.g; el = p.el;
    theta = z(1); thetaDot = z(2);
    
    thetaDoubleDot = -g/el *sin(theta);
    
    zDot = [thetaDot; thetaDoubleDot];
end
