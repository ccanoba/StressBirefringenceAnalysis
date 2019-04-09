%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the wavevector of a refracted beam, by using the
% wavevector of the incident beam, the normal vector of the surface and the
% indices of refraction of the mediums.
%
% inputs
% 
% Ki : wavevector incident beam
% NVect : Normal vector of the surface
% n1,n2 : indices of refraction of the mediums
% 
% outputs
% 
% Kr : wavevector refracted beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Kr] = SnellCalc (Ki, NVect, n1, n2)

KixN = cross(Ki,NVect);
KixN2 = KixN*n1/n2;
theta = asind(norm(KixN2));

theta = theta+180;

if theta == 180
    Kr = Ki;
else

Rotaxis = KixN/norm(KixN);

P = NVect;

Kr = (Rotaxis'*P)*Rotaxis + cosd(theta)*(P-(Rotaxis'*P)*Rotaxis)+sind(theta)*cross(Rotaxis,P);

Kr = Kr/norm(Kr);
end