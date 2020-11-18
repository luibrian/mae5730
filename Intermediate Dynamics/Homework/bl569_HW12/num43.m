%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

%43.
%Using PLOT3 (or any other way of making 3D drawings) draw a box
%with sides b,w, and h and a line through the box in direction n. 
%Then, by repeatedly using PLOT3 rotate the box smoothly rotating
%about axis n. Use this to illustrate the 120 rotation
%about the axis [1; 1; 1].

clear all;
close all;

b = 1;          %length of box base
w = 1;          %length of box width
h = 1;          %length of box height
n = [1,1,1];    %direction of line through box, angle of rotation

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


%Making the box. It will share three edges with the x,y,z axes
%Need to plot the 12 lines that make up the box
%Will have 12 edges

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
%Plotting the initial box
plot3([-2, 2], [-2, 2], [-2, 2], 'k');  %Axis of rotation
hold on;
plot3(edgeAB(1,:),edgeAB(2,:), edgeAB(3,:), 'g');
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

axis([-2 2 -2 2 -2 2]);
axis equal;
hold off;


%setting up and doing the animation
npointsT = 51;
desiredRotation =120;

%for i = 1:npointsT
for i = linspace(0,desiredRotation, npointsT)
    R = rotmat(n, i);
    %Rotate each vertex
    xNew = zeros(8,1);
    yNew = zeros(8,1);
    zNew = zeros(8,1);
    for j = 1:length(x)
        vertexOld = [x(j) y(j) z(j)]';
        newVertex = R * vertexOld;
        
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
    
    
    plot3([-2, 2], [-2, 2], [-2, 2], 'k');  %Axis of rotation
    hold on;
    plot3(edgeAB(1,:),edgeAB(2,:), edgeAB(3,:), 'g');
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

    axis([-2 2 -2 2 -2 2]);
    axis equal;
    hold off;
    pause(0.001);

end



