function [nhat, thetad] = axisAngle(R)
%AXISANGLE Takes a rotation matrix and calculates the axis and angle of
%rotation
%R is a 3x3 rotation matrix
%nhat is the axis of rotation
%thetad is the angle of rotation in degrees

Rs = (R - R')/2;
a = [-Rs(2,3); Rs(1,3); -Rs(1,2)];

nhat = a/norm(a);
stheta = norm(a);                   %sin(theta)

ctheta = (trace(R) - 1) /2;         %cos(theta)

theta = atan2(stheta, ctheta);      %in rads

thetad = theta*180/pi;

end

