function box3d(x,y,z)
%BOX3D Plots  box in 3D 
%Plots a box in 3D in which each column of x,y,z represents a vertex of
%the box. 
%x,y,z are 1x8 vectors


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



end

