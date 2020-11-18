function M_L = rhsStuffmassMatrixLagrange(I1,I2,I3,L1,L2,d1,d2,d3,m1,m2,m3,theta1,theta2,theta3)
%RHSSTUFFMASSMATRIXLAGRANGE
%    M_L = RHSSTUFFMASSMATRIXLAGRANGE(I1,I2,I3,L1,L2,D1,D2,D3,M1,M2,M3,THETA1,THETA2,THETA3)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    02-Dec-2018 01:34:47

t2 = L1.^2;
t3 = theta1-theta2;
t4 = cos(t3);
t5 = L2.*m3;
t6 = d2.*m2;
t7 = t5+t6;
t8 = L1.*t4.*t7;
t9 = theta1-theta3;
t10 = cos(t9);
t11 = L1.*d3.*m3.*t10;
t12 = theta2-theta3;
t13 = cos(t12);
t14 = L2.*d3.*m3.*t13;
M_L = reshape([I1+m2.*t2+m3.*t2+d1.^2.*m1,t8,t11,t8,I2+L2.^2.*m3+d2.^2.*m2,t14,t11,t14,I3+d3.^2.*m3],[3,3]);
