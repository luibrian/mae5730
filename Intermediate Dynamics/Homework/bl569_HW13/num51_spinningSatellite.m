%Brian Lui
%bl569
%MAE 5730- Intermediate Dynamics

clear all;
close all;

%51. Compare the three satellites:
%All three Satelllites consists of uniform-density (rho = 1) 
%objects welded together: a unit cube and spheres each with 
%radius 0.5.
%* The spheres are welded to the middles of cube faces, just 
%touching at one point.
%* In the defining configuration, the cube corners are at 
%x; y; z = 0 or 1.
%* All satellites have a massless antennas, which go from 
%(0,0,0) to (3,3,3) and the other from (0,0,0.5) to (2,2,0.5)

%Satellite A has two spheres welded to the surfaces with 
%outward normals in the +x and +y directions;
%Satellite B has a the third sphere on the +z face; and
%Satellite C has spheres on all 6 faces.

%i) Find the center of masses of the three satellites
cubeLength = 1;     %side length of cube
r = 0.5;            %radius of a sphere
rho = 1;            %Density of sphere

cm_cube = [cubeLength/2 cubeLength/2 cubeLength/2];
cm_sphere1 = [cubeLength+r, r, r]; %located in +x
cm_sphere2 = [r, cubeLength+r, r]; %located in +y
cm_sphere3 = [r, r, cubeLength+r]; %located in +z
cm_sphere4 = [-r, r, r]; %located in -x
cm_sphere5 = [r, -r, r]; %located in -y
cm_sphere6 = [r, r, -r]; %located in -z

massS = rho * (4/3*pi*r^3);     %Equals pi/6
massC = rho * cubeLength^3;      %Equals 1

%SatA CM
totMassA = 2*massS + massC;
cm_A = (massC*cm_cube + massS*cm_sphere1 + massS*cm_sphere2)./totMassA;

%SatB CM
totMassB = 3*massS + massC;
cm_B = (massC*cm_cube + massS*cm_sphere1 + massS*cm_sphere2 +...
    massS*cm_sphere3)/totMassB;

%SatC CM
totMassC = 6*massS + massC;
cm_C = (massC*cm_cube + massS*cm_sphere1 + massS*cm_sphere2 +...
    massS*cm_sphere3 + massS*cm_sphere4 +...
    massS*cm_sphere5 + massS*cm_sphere6)/totMassC;



%ii) Find the moments of inertia matrix relative to the center of mass.
%The I about the isolated object's CM
IcubeCM   = 1/6*eye(3);
IsphereCM = 2/5*massS*r^2 * eye(3);     %pi/60

%Using the parallel axis theorem to calculate
%Satellite A:
I_cube_A    = parallelAxisI(IcubeCM, massC, cm_cube, cm_A); 
I_sphere1_A = parallelAxisI(IsphereCM, massS, cm_sphere1, cm_A); 
I_sphere2_A = parallelAxisI(IsphereCM, massS, cm_sphere2, cm_A); 

I_A = I_cube_A + I_sphere1_A + I_sphere2_A;

%Satellite B:
I_cube_B    = parallelAxisI(IcubeCM, massC, cm_cube, cm_B); 
I_sphere1_B = parallelAxisI(IsphereCM, massS, cm_sphere1, cm_B); 
I_sphere2_B = parallelAxisI(IsphereCM, massS, cm_sphere2, cm_B); 
I_sphere3_B = parallelAxisI(IsphereCM, massS, cm_sphere3, cm_B); 

I_B = I_cube_B + I_sphere1_B + I_sphere2_B + I_sphere3_B;

%Satellite B:
I_cube_C    = parallelAxisI(IcubeCM, massC, cm_cube, cm_C); 
I_sphere1_C = parallelAxisI(IsphereCM, massS, cm_sphere1, cm_C); 
I_sphere2_C = parallelAxisI(IsphereCM, massS, cm_sphere2, cm_C); 
I_sphere3_C = parallelAxisI(IsphereCM, massS, cm_sphere3, cm_C); 
I_sphere4_C = parallelAxisI(IsphereCM, massS, cm_sphere4, cm_C); 
I_sphere5_C = parallelAxisI(IsphereCM, massS, cm_sphere5, cm_C); 
I_sphere6_C = parallelAxisI(IsphereCM, massS, cm_sphere6, cm_C); 

I_C = I_cube_C + I_sphere1_C + I_sphere2_C + I_sphere3_C +...
    I_sphere4_C + I_sphere5_C + I_sphere6_C; 

%I_A 
%I_B 
%I_C



%iii)
%Find the three principle moments of inertia (2 are the same for satellite A, why?)
%They are the ones aligned with the eigenvalues. 
[I_A_principleEigVect, I_A_principleEigValues] = eigs(I_A)
[I_B_principleEigVect, I_B_principleEigValues] = eigs(I_B)
[I_C_principleEigVect, I_C_principleEigValues] = eigs(I_C)



%iv) Animate these satellites tumbling in space (take care 
%that the centers of mass are stationary)
%I: rotating about each principle axis
%II: The most arbitrary motion you can find. Look at the motion
%of the antenna tips. Can you describe this motion in a simple way?
%III: If the answer to the problem above is in any way 
%simple, can you find it by pencil and paper reasoning?
 




function I0 = parallelAxisI(Icm, m, rG, rO) 
    %Given a moment of inertia matrix for an object about its CM, and
    %another point that we want to actually rotate about (in this case, the
    %CM of the entire object spheres + cube), then calculates the new
    %moment of inertia matrix about pointO
    %Icm = 3x3 matrix
    %m = mass of the object, ie of the sphere or the cube
    %rG is a 1x3 or 3x1 vector that is the object's CM
    %rO is a 1x3 or 3x1 vector that is the actual CM of the entire
    %object
    
    %making the vectors consistent
    if isrow(rO)
        rO = rO';
    end
    if isrow(rG)
        rG = rG';
    end
    
    rGrelO = rG - rO;
    
    xGrelO = rGrelO(1);
    yGrelO = rGrelO(2);
    zGrelO = rGrelO(3);
    
    offsetMat = m *...
        [(yGrelO^2 + zGrelO^2), -(xGrelO*yGrelO), -xGrelO*zGrelO;...
        -(xGrelO*yGrelO), (xGrelO^2 + zGrelO^2), -yGrelO*zGrelO;...
        -xGrelO*zGrelO, -yGrelO*zGrelO, xGrelO^2 + yGrelO^2];

    I0 = Icm + offsetMat;
end
