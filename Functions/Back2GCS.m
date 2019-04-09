%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes back all ellipsiod and beam principal directions to
% the global coordinate system.
%
% inputs
% 
% NbeamDir : principal directions index ellipsoid and beam normal modes
% Srot : rotations permformed in SetLocalCord
%
% output
% 
% A : principal directions and beam normal modes in global coordinate
% system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = Back2GCS(NbeamDir, Srot)

A = NbeamDir;

t = [0 0 -Srot(3)];

for l = 1:size(A,2)
    A(:,l) = rotateVect(A(:,l), t);
end

t = [0 -Srot(2) 0];

for l = 1:size(A,2)
    A(:,l) = rotateVect(A(:,l), t);
end

t = [-Srot(1) 0 0];

for l = 1:size(A,2)
    A(:,l) = rotateVect(A(:,l), t);
end

function [V] = rotateVect(V, t)

Rotx = @(a) [1, 0, 0 ; 0, cosd(a), -sind(a) ;0, sind(a), cosd(a)];
Roty = @(b) [cosd(b), 0, sind(b) ; 0, 1, 0 ; -sind(b), 0, cosd(b)];
Rotz = @(g) [cosd(g), -sind(g), 0 ; sind(g), cosd(g), 0 ; 0, 0, 1];

V = Rotx(-t(1))*Roty(-t(2))*Rotz(-t(3))*V;

end

end