%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots an ellipse using the parameter that describe a
% polarization ellipse
%
% Inputs :
%
% psi : Orientation angle of the ellipse
% chi : Ellipticity angle
% Makeplot : binary option for whether plot or not, used for single plots
%
% Outputs :
%
% FigID : If Makeplot is 1, return the handled of the figure, is not return
%               0
% xp, yp : vectors with the coordinates to plot the requested ellipse
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

function [FigID, xp, yp] = PlotEllip(psi,chi,Makeplot)

if nargin < 3
  Makeplot = 1;
end

% The user sees the ray comming, observation is towards -z
psi=-psi; 
t = linspace(0,2*pi,1000);
c=tan(chi);
a=sqrt(1/(1+c^2));
b=a*c;
x = a*sin(t);
y = b*cos(t);

xp=x*cos(psi)+y*sin(psi);
yp=-x*sin(psi)+y*cos(psi);
if Makeplot == 1
  plot([0,0],[-1,1],'r'),hold on,plot([-1,1],[0,0],'r')
  FigID = plot(xp,yp,'LineWidth',2.5);
  hold off
  axis ([-1 1 -1 1])
else
  FigID = 0;
end
end


 