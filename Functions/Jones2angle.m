%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calcualtes the orientation angle and ellipsticity angle
% from the Jones vector describen the polarization state of a beam
%
% Input :
%
% Jones : Jones vector of the beam
%
% Output : 
%
% psi : Orientation angle of the ellipse
% chi : Ellipticity angle
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

function [psi,chi]=Jones2angle(Jones)

Ax=sqrt(Jones(1,1)*conj(Jones(1,1)));
Ay=sqrt(Jones(2,1)*conj(Jones(2,1)));
dx=atan2(imag(Jones(1,1)),real(Jones(1,1)));
dy=atan2(imag(Jones(2,1)),real(Jones(2,1)));
d=dy-dx;
psi=1/2*atan2((2*Ax*Ay)*cos(d),(Ax^2-Ay^2));
chi=1/2*asin(((2*Ax*Ay)/(Ax^2+Ay^2))*sin(d));