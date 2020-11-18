%Brian Lui
%bl569
%MAE5730
%HW#5

clear all;
close all;

%18. 2 Masses connected by a spring

%Parameters
p.m1 = 1;                 %Mass 1
p.m2 = 1;                 %Mass 2
p.F0 = 1;                 %Magnitude of step input
p.k = 1;                  %Spring constant
p.el0 = 1;                %Rest length of the spring

%Initial condition- At these positions, stationary (velocities =0)
x1_0 = 0; x2_0 = p.el0; 
xG_0 = (1/(p.m1+p.m2))* ((p.m1*x1_0) + (p.m2*x2_0));
z0 = [x1_0; x2_0; xG_0; 0; 0; 0];    %[x1,x2,xG,v1,v2,vG]

dur = 10;
npoints = 1001;
tspan = linspace(0,dur, npoints);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part a
f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);

x1Array = zArray(:,1);
x2Array = zArray(:,2);
xGArray = zArray(:,3);
yArray = zeros(length(tArray),1);

figure(1);
axis([0, 10, 0, 45]);
hold on;
plot(tArray, x1Array,'r');
plot(tArray, x2Array,'b');
plot(tArray, xGArray,'g');
legend('m1','m2','mG');
xlabel('time')
ylabel('x position')
% 
% %%%animating 
figure(2)
for i = 1:length(tArray)
    plot(x1Array(i), yArray(i),'b*');
    hold on
    axis([-2, 30, -2, 2]);
    plot(x2Array(i), yArray(i),'r*');
    plot(xGArray(i), yArray(i), 'g*');
    xlabel('x position')
    ylabel('y position')
    hold off
    pause(0.1);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%
%part c
kArray = [1 10 100 1000 10000];
x1ErrorsArray = zeros(length(tspan),length(kArray));
x2ErrorsArray = zeros(length(tspan),length(kArray));
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
for i = 1:length(kArray)
    p.k = kArray(i);
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan, z0);
    x1Array = zArray(:,1);
    x2Array = zArray(:,2);
    xGArray = zArray(:,3);
    
    %The oscillation error in x1 is the same as x2 (rel to CM)
    x2Error = abs(x2Array - xGArray);
    x2ErrorsArray(:,i) = x2Error;
    
    figure(i+2)
    subplot(2,1,1);
    plot(tArray, x2ErrorsArray(:,i));
    xlabel('Time')
    ylabel('Error (x2-xG)')
    title(['x2-xG Error vs t, k=', num2str(kArray(i))]);
    
    subplot(2,1,2);
    plot(tArray, x1Array);
    hold on
    plot(tArray, x2Array);
    plot(tArray, xGArray);
    legend('x1','x2','xG');
    xlabel('Time')
    ylabel('Position')
    title(['Position vs Time, k=', num2str(kArray(i))]);
end

%%%%%%%%%%%%%%



function zdot = rhs(z,p)
    %Equations of motion for two masses connected by a spring
    x1 = z(1); x2 = z(2); xG = z(3); 
    x1Dot = z(4); x2Dot = z(5); xGDot = z(6);
    m1 = p.m1; m2 = p.m2; F0 = p.F0; k = p.k; el0 = p.el0;
    
    springForce = k*(x2-x1-el0);
    
    x1DoubleDot = springForce/m1;
    x2DoubleDot = F0/m2 - springForce/m2;
    xGDoubleDot = F0/(m1+m2);
    
    zdot = [x1Dot; x2Dot; xGDot; x1DoubleDot; x2DoubleDot; xGDoubleDot];
end

