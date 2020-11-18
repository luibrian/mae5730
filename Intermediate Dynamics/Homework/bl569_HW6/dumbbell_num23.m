%Brian Lui
%bl569
%MAE5730

clear all;
close all;

%Number 23- 2D dumbbell
%Two equal masses m =1 constrained by rod to be el=1 apart. There are
%non-negligible forces, except from the constraint. At t=0, they have equal
%and opposite velocites (v=1) perpendicular to the rod).

%Parameters
p.m1 = 1;                   %Mass 1
p.m2 = 1;                   %Mass 2
p.el = 1;                   %Length of the connecting rod

dur = 10;
npoints = 1001;
tspan = linspace(0, dur, npoints);

%Initial Conditions
%Can modify x1_0, y1_0, x1Dot_0, y1Dot_0, x2_0; rest is defined by
%constraint for position and velocity
x1_0 = 0.5;                 y1_0 = 0;
x1Dot_0 = 0;                y1Dot_0 = 1;
x2_0 = -0.5;                

y2_0 = sqrt(p.el^2 - (x2_0 - x1_0)^2) + y1_0;
x2Dot_0 = -x1Dot_0;         y2Dot_0 = -y1Dot_0;

z0 = [x1_0; y1_0; x1Dot_0; y1Dot_0;
    x2_0; y2_0; x2Dot_0; y2Dot_0]'; %[x1,y1,x1Dot,y1Dot, x2,y2,x2Dot,y2Dot]';

%Doing ode45 and extracting fun things from it (position of masses)
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0, options);

x1Array  = zArray(:,1);
y1Array  = zArray(:,2);
vx1Array = zArray(:,3);
vy1Array = zArray(:,4);
x2Array  = zArray(:,5);
y2Array  = zArray(:,6);
vx2Array = zArray(:,7);
vy2Array = zArray(:,8);


%Plot it!
figure(1);      %Animation of the masses in motion
title('Position in 2D of m_{1} and m_{2} connected by a rigid bar');
xlabel('x Position');
ylabel('y Position');
for i = 1:length(tArray)
    plot(x1Array(i), y1Array(i), 'g*')
    hold on;
    plot(x2Array(i), y2Array(i), 'b*')
    hold off;
    axis([-1 1 -1 1]);
    pause(0.0001)
end


%constraint check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraintError = ones(length(tArray),1)*p.el.^2 -...
    ((x2Array-x1Array).^2 + (y2Array-y1Array).^2);
figure(2)
plot(tArray, constraintError);
title('Distance between points^{2} vs Time');
xlabel('Time');
ylabel('Distance between points^{2}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Energy check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v1tSquared = vx1Array.^2 + vy1Array.^2;  %matrix
v2tSquared = vx2Array.^2 + vy2Array.^2;
v10Squared = ones(length(tArray),1) .* (x1Dot_0.^2 + y1Dot_0.^2);
v20Squared = ones(length(tArray),1) .* (x2Dot_0.^2 + y2Dot_0.^2);

KEerror = p.m1*(v1tSquared-v10Squared) + p.m2*(v2tSquared - v20Squared);

figure(3)
plot(tArray, KEerror);
title('Error in Energy Conservation vs Time');
xlabel('Time');
ylabel('Error in Energy Conservation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Angular Momentum Conservation check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angMom1 = x1Array.* vy1Array - vx1Array.*y1Array;
angMom2 = x2Array.* vy2Array - vx2Array.*y2Array;
totalAngMom = angMom1 + angMom2;
angMomInitial = ((x1_0*y1Dot_0 - x1Dot_0*y1_0) + (x2_0*y2Dot_0 - x2Dot_0*y2_0)) .* ones(length(tArray),1);
         
angMomError = angMomInitial - totalAngMom;

figure(4)
plot(tArray, angMomError);
title('Error in Angular Momemntum Conservation vs Time');
xlabel('Time');
ylabel('Error in Angular Momentum Energy Conservation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Momemntum Conservation Check
%%%%%%%%%%%%%%%%%%%%
mom = p.m1.*(sqrt(vx1Array.^2 + vy1Array.^2)) + p.m2.*(sqrt(vx2Array.^2 + vy2Array.^2));
initialMom = p.m1*(sqrt(x1Dot_0^2 + y1Dot_0^2)) + p.m2*(sqrt(x2Dot_0^2 + y2Dot_0^2));
initialMomArray = initialMom .* ones(length(tArray),1);
momError = initialMomArray - mom;

figure(5)
plot(tArray, momError);
title('Error in Linear Momemntum Conservation vs Time');
xlabel('Time');
ylabel('Error in Linear Momentum Energy Conservation');
%%%%%%%%%%%%%%%%%%%%


function zDot = rhs(z,p)
    %Getting equations of motion.
    %DAE Method- 4 from F=ma and 1 from constraint
    %Only force on masses are the constraint force imposed by the rod
    %connecting them
    x1 = z(1); y1 = z(2); x1Dot = z(3); y1Dot = z(4);
    x2 = z(5); y2 = z(6); x2Dot = z(7); y2Dot = z(8);
    
    m1 = p.m1; m2 = p.m2; el = p.el;        %parameters
    
    A = [m1, 0, 0, 0, (x1-x2)/el;
        0, m1, 0, 0, (y1-y2)/el;
        0, 0, m2, 0, (x2-x1)/el;
        0, 0, 0, m2, (y2-y1)/el;
        x2-x1, y2-y1, x1-x2, y1-y2, 0];
    b = [0;0;0;0; (x2Dot-x1Dot).^2 + (y2Dot-y1Dot).^2];
    
    %Solving Aq = b
    q = A\b;
    x1DoubleDot = q(1);
    y1DoubleDot = q(2);
    x2DoubleDot = q(3);
    y2DoubleDot = q(4);
    
    zDot = [x1Dot; y1Dot; 
        x1DoubleDot; y1DoubleDot;
        x2Dot; y2Dot; 
        x2DoubleDot; y2DoubleDot];
end
