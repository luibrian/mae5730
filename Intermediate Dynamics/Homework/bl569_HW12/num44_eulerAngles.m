%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all
close all

%Finding the net rotation of rotating about x axis, then the fixed y axis,
%then the fixed z axis, ie. Rnet = Rz(thetaz)Ry(thetay)Rx(thetax)
%Takes degrees in rads!!!!

syms thetax thetay thetaz real

thetax = thetax * pi/180;
thetay = thetay * pi/180;
thetaz = thetaz * pi/180;

Rx_thetax_fixed = [1 0 0;...
    0 cos(thetax) -sin(thetax);...
    0 sin(thetax) cos(thetax)];
Ry_thetay_fixed = [cos(thetay) 0 sin(thetay);...
    0 1 0;...
    -sin(thetay) 0 cos(thetay)];
Rz_thetaz_fixed = [cos(thetaz) -sin(thetaz) 0;...
    sin(thetaz) cos(thetaz) 0;...
    0 0 1];



Rnet_fixed = simplify(Rz_thetaz_fixed * (Ry_thetay_fixed * Rx_thetax_fixed))





%Now try doing this with the new angles
%Do this with axis angle because it's harder to figure out the components
%Rotation about ex
n1 = [0;0;1];
n1n1 = n1*n1';
skewn1 = [0, -n1(3),  n1(2);...
      n1(3),      0, -n1(1);...
     -n1(2),  n1(1),     0];
    
Rz_new = n1n1 + cos(thetaz)*(eye(3)-n1n1) + sin(thetaz)*skewn1;
%This is also equal to 
%Rz_thetaz_new = Rz_thetaz_fixed;


%Rotation about new ey
n2 = Rz_new(:,2);
n2n2 = n2*n2';
skewn2 = [0, -n2(3),  n2(2);...
      n2(3),      0, -n2(1);...
     -n2(2),  n2(1),     0];
Ry_new = (n2n2 + cos(thetay)*(eye(3)-n2n2) + sin(thetay)*skewn2)*...
    Rz_new;
Ry_new = simplify(Ry_new);
%It's also equal to 
%Ry_thetay_new = Rz_thetaz_fixed * Ry_thetay_fixed


%Rotation about new ex
n3 = Ry_new(:,1);
n3n3 = n3*n3';
skewn3 = [0, -n3(3),  n3(2);...
      n3(3),      0, -n3(1);...
     -n3(2),  n3(1),     0];
Rx_new = (n3n3 + cos(thetax)*(eye(3)-n3n3) + sin(thetax)*skewn3)*...
    Ry_new;
Rx_new = simplify(Rx_new)

%Comparing the two: If it's equal it outputs a 0 matrix and it does!
simplify(Rx_new - Rnet_fixed);


matlabFunction(Rnet_fixed, 'file', 'rotationEuler');
