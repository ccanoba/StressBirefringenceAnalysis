%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function sets index ellipsoid as global coordinates and rotates beam
% coordinate into this new coordinates system.
%
% inputs
% 
% A1, A2, A3 : principal directions index ellipsoid
% k, ke1, ke2 : beam wave vector and field directions
%
% outputs
% 
% A : rotated vectors
% Srot : angles of rotation
% Norder : order of the vector for correspondence with refractive index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Srot, Norder] = SetLocalCord(A1,A2,A3, k, ke1, ke2)

% Rotates beam coordinates to the ellipsoid coordinates

Srot = [0 0 0];
Dref = [0 0 1];     % defines the direction of propagation of light

[~,ARot] = max([Dref*A1, Dref*A2, Dref*A3]);    % find the closest component in z, that's why in stressbir all directions are oriented.
Norder(1) = ARot;   % Saves which vector goes in a particular direction

A = [A1, A2, A3, k, ke1, ke2];

t = [-atan2d(A(2, ARot),A(3, ARot)), 0, 0];

Srot(1) = -atan2d(A(2, ARot),A(3, ARot));

for l = 1:6
    A(:,l) = rotateVect(A(:,l), t);
end

t = [0, 90-acosd(A(1, ARot)), 0];

Srot(2) = 90-acosd(A(1, ARot));

for l = 1:6
    A(:,l) = rotateVect(A(:,l), t);
end

Dref = [1 0 0];
[~,ARot] = max(abs([Dref*A(:,1), Dref*A(:,2), Dref*A(:,3)]));
Norder(2) = ARot;
Norder(3) = find(~ismember([1 2 3], Norder));

t = [0, 0, atan2d(A(2, ARot),A(1, ARot))];

Srot(3) = atan2d(A(2, ARot),A(1, ARot));

for l = 1:6
    A(:,l) = rotateVect(A(:,l), t);
end

if A(2, 6)<0
    A(:,6)=A(:,6)*-1;
end

function [V] = rotateVect(V, t)

Rotx = @(a) [1, 0, 0 ; 0, cosd(a), -sind(a) ;0, sind(a), cosd(a)];
Roty = @(b) [cosd(b), 0, sind(b) ; 0, 1, 0 ; -sind(b), 0, cosd(b)];
Rotz = @(g) [cosd(g), -sind(g), 0 ; sind(g), cosd(g), 0 ; 0, 0, 1];

V = Rotx(-t(1))*Roty(-t(2))*Rotz(-t(3))*V;

end

end