%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function plots an array of polarization ellipses 
%
% Inputs :
% 
% J : It is a cell with n Jones matrices
% In : Input polarization state using Jones notation
% shift : Used to introduce a global offset in the ellipses parameters
% Efactor : Parameter used to proportionally increase the ellipticity to
%               improve visualization
% Fnum : Figure number
% px, py : Vectors with the coordinates of polarization ellipses, their
%               length equal the number of Jones matrices in J
% chiThreshold : Ellipticity angles below this value are considered
%               linearly polarization states
% ellipsize : Scales ellipse size for visualization
% arrowsize : size of the arrows in plot
% stepplot : vector with the indices of the elements of the cell J that
% will be evaluated
%
% Output
% 
% Plot with polarization ellipses.
% Red ellipses have right polarization, blue ellipses have left
% polarization, black ellipses are below threshold and considered
% quasilinear
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Jones_Ellipse_Plot(J,In,shift, Efactor, Fnum, px, py, chiThreshold, ellipsize, arrowsize, stepplot)

psishift=shift(1);
chishift=shift(2);

addpath('../Functions/ArrowHead')

for lt = 1:length(stepplot)
    
l = stepplot(lt);

Out = J{l}*In;

[psi,chi] = Jones2angle(Out);

if chi<chiThreshold
    chishiftV=0;
else
    chishiftV=chishift;
end
[~, xellip, yellip]=PlotEllip(psi+psishift,(chi+(chishiftV)*sign(chi)).*Efactor,0);
figFlag = figure(Fnum); 
hold on,
screensize = get( groot, 'Screensize' );
figFlag.Position = [screensize(3)/2-(0.8*screensize(4)/2), screensize(4)/2-(0.8*screensize(4)/2), 0.8*screensize(4), 0.8*screensize(4)];

if chi>chiThreshold
    EColor = 'r';
elseif chi<-chiThreshold
    EColor = 'b';
    xellip = fliplr(xellip);
    yellip = fliplr(yellip);
elseif abs(chi)<=chiThreshold
    EColor = 'k';
end
h(lt)=plot((xellip*ellipsize)+px(l),(yellip*ellipsize)+py(l),EColor,'LineWidth', 2);

end
axis image
line2arrow(h,'HeadLength',arrowsize,'HeadWidth',arrowsize)
end