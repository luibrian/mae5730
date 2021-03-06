function b = rhsStuffbVector(R,d,g,gamma,m,thetaDot,theta)
%RHSSTUFFBVECTOR
%    B = RHSSTUFFBVECTOR(R,D,G,GAMMA,M,THETADOT,THETA)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    07-Nov-2018 00:59:24

t2 = thetaDot.^2;
t3 = cos(theta);
b = [m.*(g.*sin(gamma)+d.*t2.*t3);-m.*(g.*cos(gamma)+d.*t2.*sin(theta));0.0;-R.*d.*m.*t2.*t3];
