%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function calculates a surface normal vectors employing the derivate
% for a pair of ortogonal directions. The normal of a surface in a
% particular point is the cross product of the lines tangential to it.
%
% Inputs
% 
% zp : matriz corresponding to the surface for normals calculation
% deltax : delta size for discrete derivate calculation
% 
% Output
%
% Normal : vector of 3 X numberofElements in the surface matrix. Every
% column correspond to a normal vector for the corresponding indented
% element of the surface matriz.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Normal] = NormalDir (zp, deltax)

[dx,dy]=gradient(zp,deltax);            % derivate of the surface in a pair of ortogonal directions

Ax = atand(dx);                                 % x coordinate tangential direction
Ay = atand(dy);                                 % y coordinate tangential direction

Vectx = [cosd(Ax(:))'; zeros(1,numel(Ax)); sind(Ax(:))'];         % Tangential vector to surface on x-z plane
Vecty = [zeros(1,numel(Ax)); cosd(Ay(:))'; sind(Ay(:))'];         % Tangential vector to surface on y-z plane

Normal = cross(Vecty,Vectx);            % Normal vector calculation

Normal = Normal./sqrt(Normal(1,:).^2+Normal(2,:).^2+Normal(3,:).^2);  % normalization of the normal vector

end