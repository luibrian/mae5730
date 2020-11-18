%Brian Lui
%bl569
%MAE5730
%HW#5

clear all;
close all;

%2 masses connected by a rigid rod

%Parameters
p.m1 = 1;                   %Mass 1
p.m2 = 1;                   %Mass 2
p.amp = 1;                  %Magnitude of sinusoidal function
p.omega = 1;                %Frequency of sinusoidal function
p.d = 1;                    %Distance between m1 and m2

dur = 10;
npoints = 1001;
tspan = linspace(0,dur, npoints);

%Initial condition: No initial velocity, two masses are d apart.
x1_0 = 0;
x2_0 = p.d;
v1_0 = 0;
v2_0 = 0;
z0 = [x1_0; x2_0; v1_0; v2_0];    

f = @(t,z) rhs(z,p,t);
[tArray, zArray] = ode45(f, tspan, z0);

%Extracting stuff from the ode45- x1, x2, calculating CM
x1Array = zArray(:,1);
x2Array = zArray(:,2);
xGArray = (p.m1+p.m2)^(-1) .* ((p.m1.*x1Array) + (p.m2.*x2Array));
xError = x2Array - xGArray;
dist = x2Array - x1Array;

%Plotting stuff
figure(1)
plot(tArray, x1Array);
hold on;
plot(tArray, x2Array);
plot(tArray, xGArray);
legend('x1','x2','xG')

figure(2)
subplot(2,1,1)
plot(tArray, xError);
title('Distance from mass 2 to Center of Mass vs Time');
xlabel('Time');
ylabel('x2 - xG');
subplot(2,1,2)
plot(tArray, dist);
title('Distance between Mass 1 and 2 vs Time');
xlabel('Time');
ylabel('x2 - x1');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trying this again but trying IC in which the masses have unequal initial
%velocities
v1_0b = 0;
v2_0b = 2;
z0b = [x1_0; x2_0; v1_0b; v2_0b];    

f = @(t,z) rhs(z,p,t);
[tArrayb, zArrayb] = ode45(f, tspan, z0b);

%Extracting stuff from the ode45- x1, x2, calculating CM
x1Arrayb = zArrayb(:,1);
x2Arrayb = zArrayb(:,2);
xGArrayb = (p.m1+p.m2)^(-1) .* ((p.m1.*x1Arrayb) + (p.m2.*x2Arrayb));
xErrorb = x2Arrayb - xGArrayb;
distb = x2Arrayb - x1Arrayb;

%Plotting stuff
figure(3)
plot(tArrayb, x1Arrayb);
hold on;
plot(tArrayb, x2Arrayb);
plot(tArrayb, xGArrayb);
legend('x1','x2','xG')



function zdot = rhs(z,p,t)
    %Equations of motion for two masses connected by a spring
    x1 = z(1); x2 = z(2); 
    x1Dot = z(3); x2Dot = z(4);
    m1=p.m1; m2=p.m2; amp=p.amp; omega=p.omega; d=p.d;
    
    F = amp*sin(omega*t);
    
    A = [m1, 0, -1; 
        0, m2, 1;
        1, -1, 0];
    b = [0;F;0];
    w = A\b;
    
    
    x1DoubleDot = w(1);
    x2DoubleDot = w(2);
    
    zdot = [x1Dot; x2Dot; x1DoubleDot; x2DoubleDot];
end


