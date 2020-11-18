function R = rotmat(n, thetad)
%ROTMAT 
%Outputs the rotation matrix about n for rotation of theta
%n is a vector
%thetad is indegrees

if isrow(n)
    n = n';
end

nhat = n/norm(n);
skew = [0,   -nhat(3),  nhat(2);...
        nhat(3),    0, -nhat(1);...
       -nhat(2), nhat(1),    0];

R1 = nhat * nhat';
R2 = cosd(thetad) * (eye(3) - R1);
R3 = sind(thetad) * skew;

R = R1 + R2 + R3;
%R = round(R);

end

