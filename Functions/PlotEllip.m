function [FigID, xp, yp] = PlotEllip(psi,chi,Makeplot)

if nargin < 3
  Makeplot = 1;
end

% chi : ellipticity --> tan(chi)=Eox/Eoy
% Eox^2+Eoy^2=1 

psi=-psi;
t = linspace(0,2*pi,1000);
c=tan(chi);
a=sqrt(1/(1+c^2));
b=a*c;
x = a*sin(t);
y = b*cos(t);

% Ex=a;
% Ey=b;
xp=x*cos(psi)+y*sin(psi);
yp=-x*sin(psi)+y*cos(psi);
if Makeplot == 1
  plot([0,0],[-1,1],'r'),hold on,plot([-1,1],[0,0],'r')
  FigID = plot(xp,yp,'LineWidth',2.5);
  % plot([-b*sin(psi),b*sin(psi)],[-b*cos(psi),b*cos(psi)],'g'),plot([-a*cos(psi),a*cos(psi)],[a*sin(psi),-a*sin(psi)],'g')
  hold off
  axis ([-1 1 -1 1])
else
  FigID = 0;
end
end


 