function b_L = rhsstuffbVectorLagrange(L1,L2,d1,d2,d3,g,m1,m2,m3,theta1,theta2,theta3,theta1dot,theta2dot,theta3dot)
%RHSSTUFFBVECTORLAGRANGE
%    B_L = RHSSTUFFBVECTORLAGRANGE(L1,L2,D1,D2,D3,G,M1,M2,M3,THETA1,THETA2,THETA3,THETA1DOT,THETA2DOT,THETA3DOT)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    02-Dec-2018 01:34:48

t2 = sin(theta1);
t3 = theta2dot.^2;
t4 = theta1-theta2;
t5 = sin(t4);
t6 = sin(theta2);
t7 = theta1dot.^2;
t8 = theta3dot.^2;
t9 = theta1-theta3;
t10 = sin(t9);
t11 = theta2-theta3;
t12 = sin(t11);
b_L = [-d1.*g.*m1.*t2-L1.*g.*m2.*t2-L1.*g.*m3.*t2-L1.*L2.*m3.*t3.*t5-L1.*d2.*m2.*t3.*t5-L1.*d3.*m3.*t8.*t10;-d2.*g.*m2.*t6-L2.*g.*m3.*t6+L1.*L2.*m3.*t5.*t7+L1.*d2.*m2.*t5.*t7-L2.*d3.*m3.*t8.*t12;d3.*m3.*(-g.*sin(theta3)+L2.*t3.*t12+L1.*t7.*t10)];
