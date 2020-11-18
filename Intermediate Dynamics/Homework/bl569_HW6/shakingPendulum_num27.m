%Brian Lui
%bl569
%MAE5730

clear all;
close all;

%Number 27- Inverted pendulum with Shaking pendulum
%Uniform stick with length el is connected to a hinge at its lower end. 
%That hind is shaking up and down with a specified acceleration a(t).

p.m = 1;                %Mass of pendulum
p.el = 10;              %Length of the pendulum
p.d = p.el/2;           %Length from base to CM of pendulum
p.g = 9.8;              %Gravity
p.omega = 500;          %Frequency of base oscillation
p.h = 0.1;              %Amplitude of base oscillation
p.I_G = 1/12 * p.m * p.el^2;                 %moment of inertia

%Time
dur = 10;
npoints = 1001;
tspan = linspace(0, dur, npoints);

%Initial Conditions
theta0 = 10* pi/180;                  %In radians
thetaDot0 = 0;
z0 = [theta0;thetaDot0];

%Doing ode45
f = @(t,z) rhs(z,p,t);
[tArray, zArray] = ode45(f, tspan, z0);

thetaArray = zArray(:,1);
xArray = -p.el .* sin(thetaArray);
yArray = p.el .* cos(thetaArray);

figure(1);
subplot(2,1,1);
plot(xArray,yArray);    
title('Position of the end of the pendulum');
xlabel('x position');
ylabel('y position');
subplot(2,1,2);
plot(tArray,thetaArray);    
title('Theta vs time');
xlabel('Time');
ylabel('Theta (rad)');

%Uncomment this to see the animation of the end of the pendulum moving in
%time
%{
figure(2);
for i = 1:length(tArray)
    plot(xArray(i), yArray(i),'*');
    hold on
    axis([-10 10 0 10]);
    title('Position of the end of the pendulum');
    xlabel('x position');
    ylabel('y position');
    pause(0.001);
end
%}


%Part c- doing it with DAE
%Doing ode45
g = @(t,z) rhsDAE(z,p,t);
[tArrayDAE, zArrayDAE] = ode45(f, tspan, z0);

thetaArrayDAE = zArrayDAE(:,1);
xArrayDAE = -p.el .* sin(thetaArrayDAE);
yArrayDAE = p.el .* cos(thetaArrayDAE);

figure(3);
subplot(2,1,1);
plot(xArrayDAE,yArrayDAE);    
title('Position of the end of the pendulum');
xlabel('x position');
ylabel('y position');
subplot(2,1,2);
plot(tArrayDAE,thetaArrayDAE);    
title('Theta vs time');
xlabel('Time');
ylabel('Theta (rad)');


%Comparing the two methods
thetaDif = thetaArray - thetaArrayDAE;
figure(4);
plot(tArray, thetaDif);
title('Difference in Theta between part b and c vs time');
xlabel('Time');
ylabel('Difference in Theta (rad)');


function zDot = rhs(z,p,t)
    %
    m=p.m; d=p.d; g=p.g; omega=p.omega; h=p.h; I_G = p.I_G;
    theta = z(1);
    thetaDot = z(2);
    
    a = -h*omega^2*cos(omega*t);  %vertical acceleration of the base due to oscillation
    
    thetaDoubleDot = (m*d*sin(theta)* (g+a)) / (I_G + m*d^2);
    zDot = [thetaDot; thetaDoubleDot];
end



function zDot = rhsDAE(z,p,t)
    %Using DAE method. The constraint is that the acceleration at O' = the
    %acceleration at A, which is the bottommost point of the pendulum.
    m=p.m; d=p.d; g=p.g; omega=p.omega; h=p.h; I_G = p.I_G;
    
    x = z(1); xDot = z(2);
    y = z(3); zDot = z(4);
    theta = z(5); thetaDot = z(6);
    
    a = -h*omega^2*cos(omega*t);  %vertical acceleration of the base due to oscillation
    
    %The DAE matrix
    A = [m, 0, 0, -1, 0;
        0, m, 0, 0, -1;
        0, 0, -I_G, d*sin(theta), d*cos(theta); 
        1, 0, d*sin(theta), 0, 0;
        0, 1, -d*cos(theta), 0, 0];
    b = [-mg; 0; 0; a-d*thetaDot^2*cos(theta); -d*sin(theta)*thetaDot];
    q = A\b;
    
    %Extracting the xDoubleDot, yDoubleDot, thetaDoubleDot
    xDoubleDot = q(1);
    yDoubleDot = q(2);
    thetaDoubleDot = q(3);
    
    zDot = [xDot; xDoubleDot; yDot; yDoubleDot; thetaDot; thetaDoubleDot];
end


