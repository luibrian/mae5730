%Brian Lui
%bl569

clear all;
close all;

%Question 32. Mass in slot on turntable

%Rigid turntable (mt ; It ) is free to rotate about a hinge at its center
%It has in it a straight frictionless slot that passes a distance 
%d from its center. A mass ms slides in the slot.

%parameters
p.ms = 1;               %Mass of moving mass
p.d  = 1;               %Distance from moving mass to center of turntable
p.It = 1;               %Moment of inertia of the turntable

v0Array = [1 2 3 2 1 -1];  %Initial velocity of moving mass
w0Array = [1 1 1 5 5 3];  %Initial rotation of the turntable

%v0 = 1;                 %initial velocity of moving mass
%w0 = 1;                 %Initial rotation of the turntable

%Legnth of simulation
dur = 100;
npoints = 1001;
tspan = linspace(0,dur,npoints);

figure(1);
hold on;
for i = 1:length(v0Array)
    v0 = v0Array(i);
    w0 = w0Array(i);
    %Initial conditions
    s0 = 0;
    sDot0 = v0;
    theta0 = 0;
    thetaDot0 = w0;
    z0 = [s0; sDot0; theta0; thetaDot0];

    %ode45
    f = @(t,z) rhs(z,p);
    [tArray, zArray] = ode45(f, tspan, z0);

    thetaArray = zArray(:,3);
    plot(tArray, thetaArray);
end
legend('v0=1,w0=1', 'v0=2,w0=1', 'v0=3,w0=1', 'v0=2,w0=5', 'v0=1,w0=5', 'v0=-1, w0=3');


function zdot = rhs(z,p)
    s = z(1); sDot = z(2); theta = z(3); thetaDot = z(4);
    ms = p.ms; d = p.d; It = p.It;
    
    thetaDDotNum = -ms*s*thetaDot*(d*thetaDot + 2*sDot);
    thetaDDotDen = ms*s^2 + It;
    thetaDoubleDot = thetaDDotNum / thetaDDotDen;
    
    sDDotNum = ms*s^3*thetaDot^2 + It*s*thetaDot^2 -...
        ms*s*thetaDot*d*(2*sDot + d*thetaDot);
    sDDotDen = ms*s^2 + It;
    sDoubleDot = sDDotNum / sDDotDen;
    
    zdot = [sDot; sDoubleDot; thetaDot; thetaDoubleDot];
end

