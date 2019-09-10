%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the optical path length for a ray traveling in a
% birefringent medium, using the average OPL. This approximation is valid
% under the assumption that there is not double refraction.
% 
% This calculation is used to determine the wavefront error due to stress
% birefringence
%
% Inputs : 
% 
% Pi, Pf : Initial and final positions of the ray, used to determine the
%            travel distance
% n1, n2 : The two refractive indices of the medium
% 
% Output :
% 
% OPL : Average optical path length for a birefringence sample
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OPL] = waveFront(Pi, Pf, n1, n2) 

L = sqrt((Pf(:,1)-Pi(:,1)).^2+(Pf(:,2)-Pi(:,2)).^2+(Pf(:,3)-Pi(:,3)).^2);   % Length of the medium
OPL = L*(n1+n2)/(2);                % Average optical path length
