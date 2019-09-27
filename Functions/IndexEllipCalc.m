%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function calculates the efective birefringence and modes of
% vibration for a beam propagating in an anisotropic media, either using
% graphical or numerical models described in literature.
%
% inputs
% 
% nx,ny,nz : components of the index ellipsoid
% S : wavevector of light
% solMethod : 1 for numerical, 2 for graphical
%
% outputs
% 
% Nbeam : indices of refraction for the components of E
% NbeamDir : directions of vibration of E
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

function [Nbeam, NbeamDir] = IndexEllipCalc(nx,ny,nz,S, solMethod)

if solMethod == 1
    [Nbeam, NbeamDir] = IndexEllipCalcNumMethod(nx,ny,nz,S);
elseif solMethod == 2
    [Nbeam, NbeamDir, OpAxis1] = IndexEllipCalcGraphMethod(nx,ny,nz,S);
else
    fprintf('Choose method 1 or 2')
end
end