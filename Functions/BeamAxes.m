%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the electric field components directions 
% for a ray with unitary wavevector k, depending on a reference 
% axis in the source coordinate system.
%
% Inputs:
%
% k : Ray unitary wavevector
% Refaxis : Reference axis in the source CS. 
%
% Outputs:
%
% ke1 : Perpendicular component of the electric field, respect to the
% reference axis
% ke2 : Parallel component of the electric field, respect to the
% reference axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors:  Camilo Cano {1*}, Pablo Zuluaga-Ramírez {2}, René Restrepo {1,3}
%   1. Applied Optics Group, Universidad EAFIT, Carrera 49 # 7 Sur-50,
%   Medellín, Colombia.
%   2. European Southern Observatory Headquarters, Karl-Schwarzschild-Str. 2, 
%   85748 Garching bei Munchen, Germany
%   3. Aerospace Optics Instrumentation Division, National Institute of Aerospace
%   Technology - INTA, Ctra de Ajalvir, Km 4, Torrejon de Ardoz, 28850 Madrid, Spain
%	* <ccanoba@eafit.edu.co>    -   2019

function [ke1, ke2] = BeamAxes(k, Refaxis)

ke1 = cross(k,repmat(Refaxis,1,size(k,2)));
ke1 = ke1./sqrt(ke1(1,:).^2+ke1(2,:).^2+ke1(3,:).^2); 
ke2 = cross(k,ke1);
ke2 = ke2./sqrt(ke2(1,:).^2+ke2(2,:).^2+ke2(3,:).^2); 

end



