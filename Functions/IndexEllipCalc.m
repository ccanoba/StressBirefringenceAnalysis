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

function [Nbeam, NbeamDir] = IndexEllipCalc(nx,ny,nz,S, solMethod)

if solMethod == 1
    [Nbeam, NbeamDir] = IndexEllipCalcNumMethod(nx,ny,nz,S);
elseif solMethod == 2
    [Nbeam, NbeamDir, OpAxis1] = IndexEllipCalcGraphMethod(nx,ny,nz,S);
else
    fprintf('Choose method 1 or 2')
end
end