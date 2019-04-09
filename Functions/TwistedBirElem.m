%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function models the Jones matrix of an anisotropic medium as a
% twisted nematic crystal.
% 
% inputs
% 
% Pi : initial position of light into the medium
% Pf : final position of light into the medium
% Bir : mean birefringence of the medium
% lambda : wavelength of the illumination
% ATi : angle of the fast axis at the input
% ATf : angle of the fast axis at the output
% 
% outputs
% 
% J : Jones matrix of the medium
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J] = TwistedBirElem(Pi,Pf,Bir,lambda, ATi, ATf)

L = sqrt((Pf(:,1)-Pi(:,1)).^2+(Pf(:,2)-Pi(:,2)).^2+(Pf(:,3)-Pi(:,3)).^2);   % Lenght of the medium

delta = (2*pi/lambda).*L.*Bir';                 % Retardance

dT = ATf-ATi;                                           % twisted angle

RotM = @(T) [cosd(T) sind(T); -sind(T) cosd(T)];        % Rotation matrix

% Equation of a twisted nematic crystal
Mret = @(Tpr, Gam) [cos(sqrt(Tpr.^2+(Gam/2)^2))+1i*Gam/2*sin(sqrt(Tpr.^2+(Gam/2)^2))/sqrt(Tpr.^2+(Gam/2)^2) Tpr*sin(sqrt(Tpr.^2+(Gam/2)^2))/sqrt(Tpr.^2+(Gam/2)^2); ...
                               -Tpr*sin(sqrt(Tpr.^2+(Gam/2)^2))/sqrt(Tpr.^2+(Gam/2)^2) cos(sqrt(Tpr.^2+(Gam/2)^2))-1i*Gam/2*sin(sqrt(Tpr.^2+(Gam/2)^2))/sqrt(Tpr.^2+(Gam/2)^2)];
                           
J = RotM(-dT)*Mret(dT*pi/180,delta); 

% Rotation of the matrix to the angle at the input ATi
J = RotM(-ATi)*J*RotM(ATi);
                           
end