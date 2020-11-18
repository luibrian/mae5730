function A = rhsStuffmassMatrixAMB(I1,I2,L1,d1,d2,m1,m2,theta1,theta2)
%RHSSTUFFMASSMATRIXAMB
%    A = RHSSTUFFMASSMATRIXAMB(I1,I2,L1,D1,D2,M1,M2,THETA1,THETA2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    20-Oct-2018 15:08:07

t2 = theta1-theta2;
t3 = cos(t2);
t4 = d2.^2;
A = reshape([-I1-L1.^2.*m2-d1.^2.*m1-L1.*d2.*m2.*t3,-L1.*d2.*m2.*t3,-I2-m2.*t4-L1.*d2.*m2.*t3,-I2-m2.*t4],[2,2]);
