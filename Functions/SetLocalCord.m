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
% beamEllipCoor : rotated vectors
% Srot : angles of rotation
% Norder : order of the vector for correspondence with refractive index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors:  Camilo Cano {1*}, Pablo Zuluaga-Ramírez {2}, René Restrepo {1,3}
%   1. Applied Optics Group, Universidad EAFIT, Carrera 49 # 7 Sur-50,
%   Medellín, Colombia.
%   2. European Southern Observatory Headquarters, Karl-Schwarzschild-Str. 2, 
%   85748 Garching bei Munchen, Germany
%   3. Aerospace Optics Instrumentation Division, National Institute of Aerospace
%   Technology - INTA, Ctra de Ajalvir, Km 4, Torrejon de Ardoz, 28850 Madrid, Spain
%	* <ccanoba@eafit.edu.co>    -   2019

function [beamEllipCoor, Srot, Norder] = SetLocalCord(A1,A2,A3, k, ke1, ke2)

% Rotates beam coordinates to the ellipsoid coordinates

Srot = [0 0 0];
Dref = [0 0 1];     % defines the direction of propagation of light, set on z axis.

[~,ARot] = max([Dref*A1, Dref*A2, Dref*A3]);    % find the closest component in z, that's why in stressbir all directions are oriented.
Norder(1) = ARot;   % Saves which vector goes in a particular direction

beamEllipCoor = [A1, A2, A3, k, ke1, ke2];

% For every direction cosine rotation to the index ellipsoid CS is calculated and defined
% by t, and then rotated using rotateVect

t = [-atan2d(beamEllipCoor(2, ARot),beamEllipCoor(3, ARot)), 0, 0];

Srot(1) = -atan2d(beamEllipCoor(2, ARot),beamEllipCoor(3, ARot));

for l = 1:6
    beamEllipCoor(:,l) = rotateVect(beamEllipCoor(:,l), t);
end

t = [0, 90-acosd(beamEllipCoor(1, ARot)), 0];

Srot(2) = 90-acosd(beamEllipCoor(1, ARot));

for l = 1:6
    beamEllipCoor(:,l) = rotateVect(beamEllipCoor(:,l), t);
end

Dref = [1 0 0];
[~,ARot] = max(abs([Dref*beamEllipCoor(:,1), Dref*beamEllipCoor(:,2), Dref*beamEllipCoor(:,3)]));
Norder(2) = ARot;
Norder(3) = find(~ismember([1 2 3], Norder));

t = [0, 0, atan2d(beamEllipCoor(2, ARot),beamEllipCoor(1, ARot))];

Srot(3) = atan2d(beamEllipCoor(2, ARot),beamEllipCoor(1, ARot));

for l = 1:6
    beamEllipCoor(:,l) = rotateVect(beamEllipCoor(:,l), t);
end

if beamEllipCoor(2, 6)<0
    beamEllipCoor(:,6)=beamEllipCoor(:,6)*-1;
end

function [V] = rotateVect(V, t)
% This function defines a rotation matrix for the transformations 
% defined by t and perform the rotation of the vector V

Rotx = @(a) [1, 0, 0 ; 0, cosd(a), -sind(a) ;0, sind(a), cosd(a)];
Roty = @(b) [cosd(b), 0, sind(b) ; 0, 1, 0 ; -sind(b), 0, cosd(b)];
Rotz = @(g) [cosd(g), -sind(g), 0 ; sind(g), cosd(g), 0 ; 0, 0, 1];

V = Rotx(-t(1))*Roty(-t(2))*Rotz(-t(3))*V;

end

end