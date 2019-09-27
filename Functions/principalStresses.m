function [principalStressDirection, principalStressValues]=PrincipalStresses(stressTensor)
%principalStresses function returns the principal stress and their
%directions using eigen function
%
% IN: stressTensor, stresses tensor to be rotated
%     
% OUT: principalStressDirection, diagonal matrix with the principal stresses
%      principalStressValues, matrix with the directions of the principal
%      stresses
%
% Authors: René Restrepo, Pablo Zuluaga y Javier Hernandez
% Grupo de optica Aplicada (Universidad EAFIT)
% Laboratorio de Instrumentación Espacial (INTA)

%% Information about principal stresses organized
[directions,values]=eig(stressTensor);

values=rot90(values);
idx = [3 2 1];
principalStressValues= values(:,idx);

principalStressDirection=directions(:,idx);
% principalStressValues=values;
% principalStressDirection=directions;
end