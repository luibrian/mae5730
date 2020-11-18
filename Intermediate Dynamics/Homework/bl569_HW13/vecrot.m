function rp = vecrot(n, thetad, r)
%VECROT Takes the vector cross product
%n is axis of rotation. It is a 3x1 vector
%thetad is the rotation angle. It is in degrees
%r is the vector to be rotated

%From the axis angle formula

if isrow(n)
    n = n';
end

nhat = n/norm(n);

r1 = nhat * dot(nhat, r);
r2 = cosd(thetad) * (r - dot(nhat,r)*nhat);
r3 = cross(nhat,r) * sind(thetad);

rp = r1 + r2 + r3;
rp = round(rp);

end

