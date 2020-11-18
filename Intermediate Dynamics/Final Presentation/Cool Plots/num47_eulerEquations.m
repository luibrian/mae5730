%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%a. 
%Draw a box that is long in the x direction and narrow in 
%the z direction so the moment of inertia in the reference
%configuration could plausibly be:
%I = [1 0 0; 0 2 0; 0 0 3];

clear all;
close all;

%Setting up I to have the desired stuff
m = 1;
b = sqrt(24/m) + 0.5;          %length of box base
w = sqrt(12/m) + 0.5;          %length of box width
h = 0.5;                       %length of box height

I0 = [w^2+h^2,       0,        0;...
           0, b^2+h^2,        0;...
           0,       0, b^2+w^2] * m/12;
I0 = round(I0);


%Making the box
vertex1 = [0, 0, 0];
vertex2 = [0, 0, h];
vertex3 = [0, w, 0];
vertex4 = [0, w, h];
vertex5 = [b, 0, 0];
vertex6 = [b, 0, h];
vertex7 = [b, w, 0];
vertex8 = [b, w, h];

%Each column represents a vertex
%    A, B, C, D, E, F, G, H
x = [0, 0, b, b, b, b, 0, 0];
y = [0, 0, 0, 0, w, w, w, w];
z = [0, h, h, 0, 0, h, h, 0];

%Edges AB, BC, CD, DA, DE, DF, FG, GH, AH, BG, CF, BE
%Bottom square- ABCD
edgeAB = [x(1), x(2); y(1), y(2); z(1), z(2)];
edgeBC = [x(2), x(3); y(2), y(3); z(2), z(3)];
edgeCD = [x(3), x(4); y(3), y(4); z(3), z(4)];
edgeAD = [x(1), x(4); y(1), y(4); z(1), z(4)];

%Top square- EFGH
edgeEF = [x(5), x(6); y(5), y(6); z(5), z(6)];
edgeFG = [x(6), x(7); y(6), y(7); z(6), z(7)];
edgeGH = [x(7), x(8); y(7), y(8); z(7), z(8)];
edgeHE = [x(5), x(8); y(5), y(8); z(5), z(8)];

%Connecting the two squares
edgeAH = [x(1), x(8); y(1), y(8); z(1), z(8)];
edgeBG = [x(2), x(7); y(2), y(7); z(2), z(7)];
edgeCF = [x(3), x(6); y(3), y(6); z(3), z(6)];
edgeDE = [x(4), x(5); y(4), y(5); z(4), z(5)];

figure(1);
plot3(edgeAB(1,:),edgeAB(2,:), edgeAB(3,:), 'g');
hold on;
plot3(edgeBC(1,:),edgeBC(2,:), edgeBC(3,:), 'b');
plot3(edgeCD(1,:),edgeCD(2,:), edgeCD(3,:), 'b');
plot3(edgeAD(1,:),edgeAD(2,:), edgeAD(3,:), 'r');
%plot top box
plot3(edgeEF(1,:),edgeEF(2,:), edgeEF(3,:), 'b');
plot3(edgeFG(1,:),edgeFG(2,:), edgeFG(3,:), 'b');
plot3(edgeGH(1,:),edgeGH(2,:), edgeGH(3,:), 'b');
plot3(edgeHE(1,:),edgeHE(2,:), edgeHE(3,:), 'b');

plot3(edgeAH(1,:),edgeAH(2,:), edgeAH(3,:), 'c');
plot3(edgeBG(1,:),edgeBG(2,:), edgeBG(3,:), 'b');
plot3(edgeCF(1,:),edgeCF(2,:), edgeCF(3,:), 'b');
plot3(edgeDE(1,:),edgeDE(2,:), edgeDE(3,:), 'b');

axis([-8 8 -8 8 -8 8]);
xlabel('x');
ylabel('y');
zlabel('z');
set(gca, 'View', [8, 35]);
hold off;



%b.    
%Simulate and animate the torque-free motion of this box using the Euler 
%equations and the initial condition w0 = [0.01, 1, 0]'. Integrate for a 
%long enough time that the motion is interesting. 
%Consider any other problems you like.

%input parameters
p.I0 = I0;      %initial inertia of object
p.M = 0;        %net moment acting on rigid object

%setting up the initial conditions
R_0 = eye(3);                      %Initally not rotating
omegaVector_0 = [0.01, 1.0, 0.0]';     %This is the initial angular velocity

%ode45 requires a column vector, so gots to convert R_0 to column
%State: R as 9x1, omegaVector as 3x1
z0 = [reshape(R_0, [9,1]); omegaVector_0];

%setting up the timescale
dur = 40;
npoints = 101;
tspan = linspace(0, dur, npoints);

%Doing ode45 to get the rotation matrix at each timestep
f = @(t,z) rhs(z,p);
[tArray, zArray] = ode45(f, tspan, z0);

%The first 9 columns make up the rotation matrix,
%the next 3 columns make up the omegaVector- don't need this for plot

%Rotating the box
for i = 1:length(tArray)
    %Reconstructing the rotation matrix
    rotAsVector = zArray(i, 1:9);
    rotAsMatrix = reshape(rotAsVector, [3,3]);
    
    %Doing the rotation
    xNew = zeros(8,1);
    yNew = zeros(8,1);
    zNew = zeros(8,1);
    for j = 1:length(x)
        vertexOld = [x(j) y(j) z(j)]';
        newVertex = rotAsMatrix * vertexOld;
        
        xNew(j) = newVertex(1);
        yNew(j) = newVertex(2);
        zNew(j) = newVertex(3);
    end
    
        %Bottom square- ABCD
    edgeAB = [xNew(1), xNew(2); yNew(1), yNew(2); zNew(1), zNew(2)];
    edgeBC = [xNew(2), xNew(3); yNew(2), yNew(3); zNew(2), zNew(3)];
    edgeCD = [xNew(3), xNew(4); yNew(3), yNew(4); zNew(3), zNew(4)];
    edgeAD = [xNew(1), xNew(4); yNew(1), yNew(4); zNew(1), zNew(4)];

    %Top square- EFGH
    edgeEF = [xNew(5), xNew(6); yNew(5), yNew(6); zNew(5), zNew(6)];
    edgeFG = [xNew(6), xNew(7); yNew(6), yNew(7); zNew(6), zNew(7)];
    edgeGH = [xNew(7), xNew(8); yNew(7), yNew(8); zNew(7), zNew(8)];
    edgeHE = [xNew(5), xNew(8); yNew(5), yNew(8); zNew(5), zNew(8)];

    %Connecting the two squares
    edgeAH = [xNew(1), xNew(8); yNew(1), yNew(8); zNew(1), zNew(8)];
    edgeBG = [xNew(2), xNew(7); yNew(2), yNew(7); zNew(2), zNew(7)];
    edgeCF = [xNew(3), xNew(6); yNew(3), yNew(6); zNew(3), zNew(6)];
    edgeDE = [xNew(4), xNew(5); yNew(4), yNew(5); zNew(4), zNew(5)];
    
    plot3(edgeAB(1,:),edgeAB(2,:), edgeAB(3,:), 'g');
    hold on;
    plot3(edgeBC(1,:),edgeBC(2,:), edgeBC(3,:), 'b');
    plot3(edgeCD(1,:),edgeCD(2,:), edgeCD(3,:), 'b');
    plot3(edgeAD(1,:),edgeAD(2,:), edgeAD(3,:), 'r');
    %plot top box
    plot3(edgeEF(1,:),edgeEF(2,:), edgeEF(3,:), 'b');
    plot3(edgeFG(1,:),edgeFG(2,:), edgeFG(3,:), 'b');
    plot3(edgeGH(1,:),edgeGH(2,:), edgeGH(3,:), 'b');
    plot3(edgeHE(1,:),edgeHE(2,:), edgeHE(3,:), 'b');

    plot3(edgeAH(1,:),edgeAH(2,:), edgeAH(3,:), 'c');
    plot3(edgeBG(1,:),edgeBG(2,:), edgeBG(3,:), 'b');
    plot3(edgeCF(1,:),edgeCF(2,:), edgeCF(3,:), 'b');
    plot3(edgeDE(1,:),edgeDE(2,:), edgeDE(3,:), 'b');

    axis([-8 8 -8 8 -8 8]);
    xlabel('x');
    ylabel('y');
    zlabel('z');

    set(gca, 'View', [8, 35]);
    hold off;
    pause(0.001);
end




function zdot = rhs(z, p)
    I0 = p.I0;
    M = p.M;
    %extracting the state
    R = reshape(z(1:9,1), [3,3]);
    omegaVector = z(10:12,1);
    
    %This is the skew symmetric version of omegaVector
    omegaTensor = [0, -omegaVector(3),  omegaVector(2);...
      omegaVector(3),               0, -omegaVector(1);...
     -omegaVector(2),  omegaVector(1),              0];
    
    I = R * I0 * R';
    omegaVectorDot = inv(I) * (eye(3,1)*M - cross(omegaVector, I*omegaVector));
    Rdot = omegaTensor * R;
    zdot = [reshape(Rdot, [9,1]); omegaVectorDot];
end