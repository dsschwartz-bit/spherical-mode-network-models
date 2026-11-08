function [Fprime,Xprime,Yprime,Zprime] = RotateVectorField(Axis,angle,X,Y,Z,F)
%ROTATEVECTORFIELD  Rotate a sampled 3-D vector field and its sample points.
%   [Fprime,Xprime,Yprime,Zprime] = RotateVectorField(Axis,angle,X,Y,Z,F)
%   rotates by angle (radians) about the 3-element direction Axis (see
%   rotationmat3D). X,Y,Z are 3-D meshgrid arrays of the sample points and
%   F is (size(X)) x 3 holding the Cartesian vector components. Each point
%   and each field vector is rotated; outputs have the same sizes as the
%   inputs.

Axis=reshape(Axis,3,1);
Rot=rotationmat3D(angle,Axis);

if any(size(X)~=size(Y))
    error("X,Y,Z need to be a meshgrid")
elseif any(size(Y)~=size(Z))
    error("X,Y,Z need to be a meshgrid")
end

dimx=size(F,1);
dimy=size(F,2);
dimz=size(F,3);
Xprime=zeros(dimx,dimy,dimz);Yprime=Xprime;Zprime=Xprime;
Fprime=zeros(dimx,dimy,dimz,3);
for n=1:dimx
    for m=1:dimy
        for o=1:dimz
            rprime=Rot*[X(n,m,o);Y(n,m,o);Z(n,m,o)];
            Xprime(n,m,o)=rprime(1);Yprime(n,m,o)=rprime(2);Zprime(n,m,o)=rprime(3);
            Fprime(n,m,o,:) = Rot*squeeze(F(n,m,o,:));
        end
    end
end

end