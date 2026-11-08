function R= rotationmat3D(angle,Axis)
%ROTATIONMAT3D  3x3 matrix for a rotation about an arbitrary axis.
%   R = rotationmat3D(angle,Axis) returns R such that R*x rotates the
%   column vector x by angle (radians, right-hand rule) about the line
%   through the origin and the point Axis (3-element, non-zero).
%
%   Examples:
%     r = rand(4,1);   % rotate and rotate back -> identity
%     rotationmat3D(r(1),r(2:4)) * rotationmat3D(-r(1),r(2:4))
%     Rz45 = rotationmat3D(pi/4,[0 0 1]);   % 45 deg about z
%
%   Bileschi 2009

L = norm(Axis);
if (L < eps)
   error('axis direction must be non-zero vector');
end
Axis = Axis / L;

u = Axis(1);
v = Axis(2);
w = Axis(3);
u2 = u^2;
v2 = v^2;
w2 = w^2;
c = cos(angle);
s = sin(angle);
R = zeros(3);
R(1,1) =  u2 + (v2 + w2)*c;
R(1,2) = u*v*(1-c) - w*s;
R(1,3) = u*w*(1-c) + v*s;
R(2,1) = u*v*(1-c) + w*s;
R(2,2) = v2 + (u2+w2)*c;
R(2,3) = v*w*(1-c) - u*s;
R(3,1) = u*w*(1-c) - v*s;
R(3,2) = v*w*(1-c)+u*s;
R(3,3) = w2 + (u2+v2)*c;