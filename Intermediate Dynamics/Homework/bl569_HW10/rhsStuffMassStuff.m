function A = rhsStuffMassStuff(I,R,d,m,theta)
%RHSSTUFFMASSSTUFF
%    A = RHSSTUFFMASSSTUFF(I,R,D,M,THETA)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    07-Nov-2018 00:59:24

t2 = sin(theta);
A = reshape([m,0.0,1.0,-m.*(R-d.*t2),-d.*m.*t2,-d.*m.*cos(theta),-R,I-d.^2.*m+R.*d.*m.*t2,1.0,0.0,0.0,-R+d.*t2,0.0,-1.0,0.0,0.0],[4,4]);