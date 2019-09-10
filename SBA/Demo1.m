%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program calculates changes in the state of polarization of light 
% after interacting with a loaded optical component. The fundamental 
% concept used in the calculations is the photoelastic law, which allows
% the prediction of changes in the birefringence of a component from 
% previous knowledge of the stress distribution it presents.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Samples')
addpath(genpath('../Functions'))

StressesDataFileName = {'SqWindow4.txt', 'lensConv.txt',...
                                        'lensConv4.txt','lensConvPress.txt', 'Lens25.txt','lensConvTest.txt','MeshT_12.txt'};

Data=load(StressesDataFileName{2});

verbosity = 1;

%% Parameters definition

x = Data(:,2)+Data(:,11); y = Data(:,3)+Data(:,12); z = Data(:,4)+Data(:,13); % Deformed nodes coordinates definition
thetaref = 0;                                                                                   % Orientation of the reference plane of polarization respect to x
k11 = -0.5e-12;                                                         % Stress optical coeficient 1
k12 = -3.3e-12;                                                         % Stress optical coeficient 1 3.3
OSC = [k11, k12];
lambda = 532e-9;                                                    % Light wavelenght
n0 = 1.5;                                               % refractive index without load
Nsx = 101;                                                                   % number of divisions in grid for interpolation
solMethod = 1;                                      % Method for solution. 1 numerical, 2 graphical.
considerDiattenuation = 1;

illumParam = sourceDefinition(2, 12e-3, 25, 36*2, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);   
% illumParam = sourceDefinition(1, 12e-3, 25, 36*2, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);   
% illumParam = sourceDefinition(6, 12e-3, 24e-3, 24e-3, 26, 26, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);  

[RetardanceJonesMatrix, Layer, beamLoc, birefringence, axesRot, WF, diattData] = StressBirRayTracing(Data, n0, OSC, lambda, illumParam, Nsx, thetaref, solMethod, considerDiattenuation, verbosity);

JonesMatrix = cell(1, length(RetardanceJonesMatrix));
for l=1:length(RetardanceJonesMatrix)
    JonesMatrix{l} = diattData{1}{l,2}*RetardanceJonesMatrix{l}*diattData{1}{l,1};
end
%%
if verbosity == 1
    Intheta = 0;
    In = [cosd(Intheta) sind(Intheta)]'; In=In/norm(In);
    shift = [0 0];
    Efactor = 1;
    OutLayer = length(Layer) ;
    l = 250;
    chiThreshold = pi/(2^7);
    ellipsize = 3;
    arrowsize = 5;
    stepplotR = 1:3:illumParam.NRI;
    stepplotT = 1:3:illumParam.NTI;
    stepplot = (reshape(repmat(stepplotR,length(stepplotT),1),length(stepplotR)*length(stepplotT),1)'-1)*illumParam.NTI+repmat(stepplotT,1,length(stepplotR));
    PosFactor = 5000;
    NumF = 26;
    
    Jones_Ellipse_Plot(JonesMatrix,In,shift, Efactor, NumF, beamLoc{OutLayer}(:,1)*PosFactor, beamLoc{OutLayer}(:,2)*PosFactor, chiThreshold, ellipsize, arrowsize, stepplot)
    hold off
    title('Polarization map')
end
%%

xs = max([abs(min(Data(:,2))), max(Data(:,2)), abs(min(Data(:,3))), max(Data(:,3))]);                 % x size of grid for interpolation
xs = linspace(-xs, xs, Nsx);                                           % x axis definition
[xs, ys] = meshgrid(xs, xs);                                           % grid creation

if verbosity == 1
    birefringenceMap = birefringence{1};
    birefringenceMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),birefringenceMap,xs,ys);
    figure,imagesc(birefringenceMap),colorbar, colormap ('plasma'), title('Front birefringence map'), colorbar
    birefringenceMap = birefringence{length(Layer)};
    birefringenceMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),birefringenceMap,xs,ys);
    figure,imagesc(birefringenceMap*-1),colorbar, colormap ('plasma'), title('Rear birefringence map'), colorbar
    axesRotMap = axesRot{1};
    axesRotMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),axesRotMap,xs,ys);
    figure,imagesc(axesRotMap),colorbar, colormap(flipud(cmap('c3','shift',0.25))), title('Front fast axis orientation map'), colorbar
    axesRotMap = axesRot{length(Layer)};
    axesRotMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),axesRotMap,xs,ys);
    figure,imagesc(axesRotMap),colorbar, colormap(flipud(cmap('c3','shift',0.25))), title('Rear fast axis orientation map'), colorbar
end

JM = cell2mat(RetardanceJonesMatrix);

retardance = 2*acos(real(JM(1,1:2:end)));
thetaEnd1 = imag(JM(2,1:2:end))./sin(retardance/2);
thetaEnd2 = imag(JM(1,1:2:end))./sin(retardance/2);
thetaEnd =atan2d(thetaEnd1, -thetaEnd2)/2;
retardanceMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),retardance',xs,ys);
thetaEndMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),thetaEnd,xs,ys);

if verbosity == 1
    figure,imagesc(WF{2,end}), colormap(colorcet('cbl1')), title('Stress birefringence wavefront error'), colorbar
    figure,imagesc(retardanceMap), colormap ('plasma'), title('Retardance'), colorbar
    figure,imagesc(thetaEndMap), colormap(flipud(cmap('c3','shift',0.25))), title('Effective fast axis orientation'), colorbar
    if considerDiattenuation == 1
        figure,imagesc(diattData{2}), colormap(flipud(cmap('c3','shift',0.25))),title('S-component orientation angle map'), colorbar
        figure,imagesc(diattData{3}), colormap('inferno'),title('Diattenuation map'), colorbar
    end
end

save('../Output/demo1Output','beamLoc','JonesMatrix','birefringence','axesRot')