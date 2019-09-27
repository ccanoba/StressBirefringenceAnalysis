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

% Authors:  Camilo Cano {1*}, Pablo Zuluaga-Ramírez {2}, René Restrepo {1,3}
%   1. Applied Optics Group, Universidad EAFIT, Carrera 49 # 7 Sur-50,
%   Medellín, Colombia.
%   2. European Southern Observatory Headquarters, Karl-Schwarzschild-Str. 2, 
%   85748 Garching bei Munchen, Germany
%   3. Aerospace Optics Instrumentation Division, National Institute of Aerospace
%   Technology - INTA, Ctra de Ajalvir, Km 4, Torrejon de Ardoz, 28850 Madrid, Spain
%	* <ccanoba@eafit.edu.co>    -   2019

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