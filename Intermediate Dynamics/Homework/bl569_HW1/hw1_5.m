%Brian Lui
%Bl569
%MAE5780
%9/5/18

%5. Simple animation. Draw a picture of some object (a face, a house, 
%whatever), and make it move around on the screen in a smooth and 
%interesting way. No distortions. Just motions.

%Creating the image (smiley face!)
points = 101;
radiusLarge = 3;
radiusSmall = .5;

A1 = circle(0,0,radiusLarge,points);
A2 = circle(-1,1,radiusSmall,points);
A3 = circle(1,1,radiusSmall,points);
A4 =semicircle(0,0,2,points,false);
liftpen = [0, inf];
smiley = vertcat(A1,liftpen, A2,liftpen,A3,liftpen,A4,liftpen)';


%define fancy functions
p_5.b = 0.1;

dur5 = 500; 
npoints_5 = 1001;
tspan_5 = linspace(0,dur5,npoints_5);
z0_5 = [0;5;3];

f_5 = @(t,z) rhs5(z, p_5);
[tArray_5, zArray_5] = ode45(f_5, tspan_5, z0_5);

xArray_5 = zArray_5(:,1);
yArray_5 = zArray_5(:,2);
zArray_5 = zArray_5(:,3);

figure(4)
% Moving the image around
% Uses the Thomas' cyclically symmetric attractor
for n = 1:npoints
    offset = [xArray_5(n); yArray_5(n)]; %add to every element on image    
    newSmiley = offset + smiley; 

    xCoord = newSmiley(1,:); 
    yCoord = newSmiley(2,:);
    plot(xCoord, yCoord);
    axis([-10 10 -10 10]);
    pause(0.0001)
end


function A = circle(x,y,r,np)
    %Returns matrix A containing x and y coordinates for a circle with 
    %center located at (x,y) and having radius r using np points
    theta = linspace(0,360, np)';
    xPlot = r*cosd(theta) + x;
    yPlot = r*sind(theta) + y;
    A = [xPlot, yPlot];
end


function A = semicircle(x,y,r,np,up)
    %Returns matrix A containing x and y coordinates for a circle with 
    %center located at (x,y) and having radius r using np points. If 
    %up = true, then semicircle  from 0-180degree, otherwise a semicircle 
    %from 180-360 degrees.
    if up == true
        theta = linspace(0,180,np)';
    else
        theta = linspace(180,360,np)';
    end
    xPlot = r*cosd(theta) + x;
    yPlot = r*sind(theta) + y;
    A = [xPlot yPlot];
end


function rDot = rhs5(r,p)
    %Thomas' cyclically symmetric attractor
    %https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor
    x = r(1);
    y = r(2);
    z = r(3);
    b = p.b;
    
    xDot = sin(y) - b*x;
    yDot = sin(z) - b*y;
    zDot = sin(x) - b*z;
    rDot = [xDot; yDot; zDot];
    
end