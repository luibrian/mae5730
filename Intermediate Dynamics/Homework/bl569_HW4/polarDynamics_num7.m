%Brian Lui
%MAE 5730
%HW4
%9/26/18

close all;
clear all;

%Simplest dynamics with Polar Coordinates
%Assume a particle is on a plane with no force on it. 
%So, you know it moves at constant speed in a constant direction (a = 0)

%b) Solve differential equations numerically for various initial conditions.
%If rDot and thetaDot both equal 0, it's just a point
%r cannot equal 0 initiall
IC = [1 0 1 0; 
    1 0 1 1;
    1 1 1 0;
    3 1 0 1;
    1 2 0 2;
    1 5 1 1;
    1 1 3 4;
    3 1 2 1;
    2 4 3 1];

dur = 10; 
npoints = 1001;
tspan = linspace(0,dur, npoints);

tolerances = [1E-2, 1E-4, 1E-6, 1E-8, 1E-10, 1E-12];

%each column of this represents a tolerance level, each row is a different
%IC

errors = zeros(length(IC(:,1)),length(tolerances));


for a = 1:length(tolerances)
    options = odeset('RelTol',tolerances(a),'AbsTol',tolerances(a));
    figure(a);  %each figure produces is a different tolerance
    hold on
    for i = 1:length(IC(:,1))
        z0 = IC(i,:)';

        f = @(t,z) rhs(z);
        [tArray, zArray] = ode45(f, tspan, z0,options);
        rArray = zArray(:,1);
        thetaArray = zArray(:,2);

        xArray = rArray.*cos(thetaArray);
        yArray = rArray.*sin(thetaArray);


        plot(xArray, yArray);
        errors(i,a) = isStraight(xArray, yArray);
    end
    hold off;
end

errors


%c) Plot the solution and check that the motion is a straight line at constant speed.

%d) Using your numerical result, pick a way to measure how straight the 
%path is, and see how straight a line your polar coordinate solution gives. 
%You should define a quantitative measure of straightness, and 
%then measure it with your solution.

%e) Is the path more straight when you redefne the numerical tolerances.

function zDot = rhs(z)
    %defining equations of motion for this case. 
    r = z(1); theta = z(2); rDot = z(3); thetaDot = z(4);
    
    rDoubleDot = r*thetaDot^2;
    thetaDoubleDot = -2*rDot*thetaDot/r;
    
    zDot = [rDot; thetaDot; rDoubleDot; thetaDoubleDot];
end

function error = isStraight(xArray,yArray)
    %takes in two vectors corresponding to the x and y coords
    %returns the difference between the xArray and the yArray 
    %if the error is 0, then it is perfectly straight
    dY = yArray(2)-yArray(1);
    dX = xArray(2)-xArray(1);
    
    deltaY = yArray(end)-yArray(1);
    deltaX = xArray(end)-xArray(1);
    
    dydx = dY/dX;
    deltaYdeltaX = deltaY/deltaX;
    
    error = abs(deltaYdeltaX - dydx);
end