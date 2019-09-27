%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes back all ellipsiod and beam principal directions to
% the global coordinate system.
%
% inputs
% 
% NbeamDir : principal directions index ellipsoid and beam normal modes
% Srot : rotations performed in SetLocalCord
%
% output
% 
% globalDCoor : principal directions and beam normal modes in global coordinate
% system
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

function [globalDCoor] = Back2GCS(NbeamDir, Srot)

% Independently rotates every axes, in opposite direction to SetLocalCord
globalDCoor = NbeamDir;

t = [0 0 -Srot(3)];

for l = 1:size(globalDCoor,2)
    globalDCoor(:,l) = rotateVect(globalDCoor(:,l), t);
end

t = [0 -Srot(2) 0];

for l = 1:size(globalDCoor,2)
    globalDCoor(:,l) = rotateVect(globalDCoor(:,l), t);
end

t = [-Srot(1) 0 0];

for l = 1:size(globalDCoor,2)
    globalDCoor(:,l) = rotateVect(globalDCoor(:,l), t);
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