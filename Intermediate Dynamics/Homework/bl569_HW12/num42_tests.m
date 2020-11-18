%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%Testing function vecrot and rotmat. Parts a,b, and c

ex = [1,0,0]';
ey = [0,1,0]';
ez = [0,0,1]';
a  = [1,1,1]';
b  = [2,2,2]';

%Checking vecrot
ntestVecrot = 15;
testVecrot = zeros(ntestVecrot, 1);

%Rotation of 90 degrees about x axis
testVecrot(1) = isequal(vecrot(ex, 90, ex), ex);
testVecrot(2) = isequal(vecrot(ex, 90, ey), ez);
testVecrot(3) = isequal(vecrot(ex, 90, ez), -ey);

%Rotation of 90 degrees about y axis
testVecrot(4) = isequal(vecrot(ey, 90, ex), -ez);
testVecrot(5) = isequal(vecrot(ey, 90, ey), ey);
testVecrot(6) = isequal(vecrot(ey, 90, ez), ex);

%Rotation of 90 degrees about z axis
testVecrot(7) = isequal(vecrot(ez, 90, ex), ey);
testVecrot(8) = isequal(vecrot(ez, 90, ey), -ex);
testVecrot(9) = isequal(vecrot(ez, 90, ez), ez);

%Rotation of 90 degrees about a axis
testVecrot(10) = isequal(vecrot(a, 120, ex), ey);
testVecrot(11) = isequal(vecrot(a, 120, ey), ez);
testVecrot(12) = isequal(vecrot(a, 120, ez), ex);

%Rotation of 90 degrees about a axis
testVecrot(13) = isequal(vecrot(b, 120, ex), ey);
testVecrot(14) = isequal(vecrot(b, 120, ey), ez);
testVecrot(15) = isequal(vecrot(b, 120, ez), ex);



%Checking rotmat
ntestAxisAngle = 12;
testRotmat = zeros(ntestAxisAngle, 1);

%Rotation of 90 degrees about x axis
Rx = rotmat(ex,90);
testRotmat(1) = isequal(Rx*ex, ex);
testRotmat(2) = isequal(Rx*ey, ez);
testRotmat(3) = isequal(Rx*ez, -ey);

%Rotation of 90 degrees about y axis
Ry = rotmat(ey,90);
testRotmat(4) = isequal(Ry*ex, -ez);
testRotmat(5) = isequal(Ry*ey, ey);
testRotmat(6) = isequal(Ry*ez, ex);

%Rotation of 90 degrees about z axis
Rz = rotmat(ez,90);
testRotmat(7) = isequal(Rz*ex, ey);
testRotmat(8) = isequal(Rz*ey, -ex);
testRotmat(9) = isequal(Rz*ez, ez);

%Rotation of 90 degrees about z axis
Ra = round(rotmat(a,120));
testRotmat(10) = isequal(Ra*ex, ey);
testRotmat(11) = isequal(Ra*ey, ez);
testRotmat(12) = isequal(Ra*ez, ex);

%Rotation of 90 degrees about z axis
Rb = round(rotmat(b,120));
testRotmat(13) = isequal(Rb*ex, ey);
testRotmat(14) = isequal(Rb*ey, ez);
testRotmat(15) = isequal(Rb*ez, ex);


testbc = [testVecrot, testRotmat]


%Part d
%Write a Matlab script that takes rotation matrix as input and calculates
%the axis and angle of rotation. Check that this works forward and back
%with your code from above, both ways, with some odd examples

%Checking rotmat
ntestAxisAngle = 10;
testAxisAngle = zeros(ntestAxisAngle, 1);

[nRx,angleRx] = axisAngle(Rx);
[nRy,angleRy] = axisAngle(Ry);
[nRz,angleRz] = axisAngle(Rz);
[nRa,angleRa] = axisAngle(Ra);
[nRb,angleRb] = axisAngle(Rb);

angleRa = round(angleRa);
angleRb = round(angleRb);

testAxisAngle(1) = isequal(nRx, ex);
testAxisAngle(2) = isequal(angleRx, 90);
testAxisAngle(3) = isequal(nRy, ey);
testAxisAngle(4) = isequal(angleRy, 90);
testAxisAngle(5) = isequal(nRz, ez);
testAxisAngle(6) = isequal(angleRz, 90);
testAxisAngle(7) = isequal(nRa, a/norm(a));
testAxisAngle(8) = isequal(angleRa, 120);
testAxisAngle(9) = isequal(nRb, b/norm(b));
testAxisAngle(10) = isequal(angleRb, 120);

testAxisAngle

