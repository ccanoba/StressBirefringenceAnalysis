%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program calculates changes in the state of polarization of light 
% after interacting with a loaded optical component. The fundamental 
% concept used in the calculations is the photoelastic law, which allows
% the prediction of changes in the birefringence of a component from 
% previous knowledge of the stress distribution it presents.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path to functions
addpath(genpath('..'))

stressesDataFileName = {'lensConv.txt', 'lensConv4.txt'};

% Load data from FEM analysis
data=load(stressesDataFileName{1});

% 1 for plotting results
verbosity = 1;

%% Parameters definition

x = data(:,2)+data(:,11); y = data(:,3)+data(:,12); z = data(:,4)+data(:,13); % Deformed nodes coordinates definition
thetaref = 0;                                                                                   % Orientation of the reference plane of polarization respect to x
k11 = -0.5e-12;                                                         % Stress optical coeficient 1
k12 = -3.3e-12;                                                         % Stress optical coeficient 1
OSC = [k11, k12];
lambda = 532e-9;                                                    % Light wavelenght
n0 = 1.5;                                               % refractive index without load
Nsx = 101;                                                                   % number of divisions in grid for interpolation
solMethod = 1;                                      % Method for solution. 1 numerical, 2 graphical.
considerDiattenuation = 1;                  % 1 to calculate diattenuation on front and rear surfaces

% Defines the model source, help SourceDefinition to identify parameters

% Collimated spherical discretization
% illumParam = SourceDefinition(2, 12e-3, 25, 36*2, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);   
% Diverging spherical discretization
% illumParam = SourceDefinition(1, 12e-3, 25, 36*2, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]); 
% Collimated square discretization. Circular shape
illumParam = SourceDefinition(6, 12e-3, 24e-3, 24e-3, 26, 26, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);  

% Ray tracing in stress birefringent component
[retardanceJonesMatrix, Layer, beamLoc, kBeam, birefringence, axesRot, WF, diattData] = StressBirRayTracing(data, n0, OSC, lambda, illumParam, Nsx, thetaref, solMethod, considerDiattenuation, verbosity);

% If diattenuation is calculated, it is multiplied by the retardance
% considering that front surfaces effect is before retardance and rear
% surface effect act after retardance
if considerDiattenuation == 1
    JonesMatrix = cell(1, length(retardanceJonesMatrix));
    for l=1:length(retardanceJonesMatrix)
        JonesMatrix{l} = diattData{1}{l,2}*retardanceJonesMatrix{l}*diattData{1}{l,1};
    end
else
    JonesMatrix = retardanceJonesMatrix;
end
%% Plots polarization map

if verbosity == 1 && ismember(illumParam.illumCase, [1,2])
    Intheta = 0;        % Angle of orientation for linear polarizarion state
    In = [cosd(Intheta) sind(Intheta)]'; In=In/norm(In);    % Define linear polarization state using Intheta
    shift = [0 0];      % Used to introduce a global offset in the ellipses parameters
    eFactor = 1;        % Parameter used to proportionally increase the ellipticity
    outLayer = length(Layer) ;
    chiThreshold = pi/(2^7);    % Ellipticity angles below this value are considered linearly polarization states
    ellipsize = 3;      % Scales ellipse size for visualization
    arrowsize = 5;
    stepplotR = 1:3:illumParam.NRI;
    stepplotT = 1:3:illumParam.NTI;
    % vector with the indices of the rays that will be plotted
    stepplot = (reshape(repmat(stepplotR,length(stepplotT),1),length(stepplotR)*length(stepplotT),1)'-1)*illumParam.NTI+repmat(stepplotT,1,length(stepplotR));
    PosFactor = 5000;   % Used to scale the coordinate of polarization ellipses
    Fnum = 26;      % Figure number
    
    Jones_Ellipse_Plot(JonesMatrix,In,shift, eFactor, Fnum, beamLoc{outLayer}(:,1)*PosFactor, beamLoc{outLayer}(:,2)*PosFactor, chiThreshold, ellipsize, arrowsize, stepplot)
    hold off
    title('Polarization map')
end
%% Visualize results

% Define space to interpolate results
xs = max([abs(min(data(:,2))), max(data(:,2)), abs(min(data(:,3))), max(data(:,3))]);                 % x size of grid for interpolation
xs = linspace(-xs, xs, Nsx);                                           % x axis definition
[xs, ys] = meshgrid(xs, xs);                                           % grid creation

if verbosity == 1
    birefringenceMap = birefringence{1};
    birefringenceMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),birefringenceMap,xs,ys);
    figure,imagesc(rot90(birefringenceMap*-1,2)),colorbar, colormap ('plasma'), title('Front birefringence map'), colorbar
    birefringenceMap = birefringence{length(Layer)};
    birefringenceMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),birefringenceMap,xs,ys);
    figure,imagesc(rot90(birefringenceMap*-1,2)),colorbar, colormap ('plasma'), title('Rear birefringence map'), colorbar
    axesRotMap = axesRot{1};
    axesRotMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),axesRotMap,xs,ys);
    figure,imagesc(rot90(axesRotMap,2)),colorbar, colormap(flipud(cmap('c3','shift',0.25))), title('Front fast axis orientation map'), colorbar
    axesRotMap = axesRot{length(Layer)};
    axesRotMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),axesRotMap,xs,ys);
    figure,imagesc(rot90(axesRotMap,2)),colorbar, colormap(flipud(cmap('c3','shift',0.25))), title('Rear fast axis orientation map'), colorbar
end

% Calculates retardance and orientation maps
JM = cell2mat(retardanceJonesMatrix);
retardance = 2*acos(real(JM(1,1:2:end)));
thetaEnd1 = imag(JM(2,1:2:end))./sin(retardance/2);
thetaEnd2 = imag(JM(1,1:2:end))./sin(retardance/2);
thetaEnd =-atan2d(thetaEnd1, -thetaEnd2)/2;
retardanceMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),retardance',xs,ys);
thetaEndMap = griddata(beamLoc{length(Layer)+1}(:,1),beamLoc{length(Layer)+1}(:,2),thetaEnd,xs,ys);

if verbosity == 1
    figure,imagesc(rot90(WF{2,end},2)), colormap(colorcet('cbl1')), title('Stress birefringence wavefront error'), colorbar
    figure,imagesc(rot90(retardanceMap,2)), colormap ('plasma'), title('Retardance'), colorbar
    figure,imagesc(rot90(thetaEndMap,2)), colormap(flipud(cmap('c3','shift',0.25))), title('Effective fast axis orientation'), colorbar
    if considerDiattenuation == 1
        figure,imagesc(rot90(diattData{2},2)), colormap(flipud(cmap('c3','shift',0.25))),title('S-component orientation angle map'), colorbar
        figure,imagesc(rot90(diattData{3},2)), colormap('inferno'),title('Diattenuation map'), colorbar
    end
end

% Save outputs for further analysis
save('../Output/demo1Output','beamLoc','kBeam','JonesMatrix','birefringence','axesRot')