%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This program calculates changes in the state of polarization of light 
% after interacting with a loaded optical component. The fundamental 
% concept used in the calculations is the photoelastic law, which allows
% the prediction of changes in the birefringence of a component from 
% previous knowledge of the stress distribution it presents.
% 

% How to run

% Add paths with examples and the functions this script employs

% If linux or windows
addpath('../Samples')
addpath('../Functions')

StressesDataFileName = {'Press_430kPa_z10.txt', 'NewModel2.txt', 'SqWindow4.txt'...
                                        'BK7_Test3_300_2.txt', 'Data1.txt', 'Lens1.txt', 'Lens2.txt','lensConv.txt','lensConv4.txt','lensConvPress.txt'};

Data=load(StressesDataFileName{3});

verbosity = 1;

%% Parameters definition

xu = Data(:,2); yu = Data(:,3); zu = Data(:,4);                                              % Undeformed nodes coordinates definition
x = Data(:,2)+Data(:,11); y = Data(:,3)+Data(:,12); z = Data(:,4)+Data(:,13); % Deformed nodes coordinates definition

% x = Data(:,2);
% y = Data(:,3);
% z = Data(:,4);

% Param = struct('RI',12e-3, 'NRI',15,'NTI',36,'P0',[0 0 -200e-3], 'Poslim', [min(x), max(x), min(y), max(y), min(z), max(z)],...
%                 'xside', 24e-3,'yside', 24e-3,'Nx', 50, 'Ny', 50);
% illumCase = 5;    % Casos 1,3,5        
            
            
% Param = struct('RI',12e-3, 'NRI',15,'NTI',36,'P0',[0 0 -200e-3], 'Poslim', [min(x), max(x), min(y), max(y), min(z), max(z)],...
%                 'xside', 49e-3,'yside', 49e-3,'Nx', 50, 'Ny', 50);
% illumCase = 8;     % Casos 2, 4

% Lens example
% illumParam = sourceDefinition(2, 12e-3, 15, 36, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);                             % 1,2

%Square example
illumParam = sourceDefinition(5, 49e-3, 49e-3, 50, 50, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);                  % 3,5

% illumParam = sourceDefinition(4, 12e-3, 49e-3, 49e-3, 50, 50, [0 0 -200e-3], [min(x), max(x), min(y), max(y), min(z), max(z)]);       % 4,6

%% Discretize nodes in the model into layers

fprintf('Nodes cloud discretization\n')

LayerFound = [];
cont = 1;
haveNodes = true;

while haveNodes==true
    Layer{cont}=FindNodeLayer(Data, LayerFound, illumParam);
    LayerFound = [LayerFound; Layer{cont}];
    cont = cont+1;
    if length(LayerFound)==length(Data)
        haveNodes = false;
    end
end

clear LayerFound cont haveNodes

%% Plot layers
% close all

if verbosity == 1
    colores = ['r','b','g'];
    for i=1:length(Layer)
        hold on,plot3(z(Layer{i}),y(Layer{i}),x(Layer{i}),[colores(mod(i,3)+1) '*'])
        xlabel('z'), ylabel('y'), zlabel('x')
        view([0.5 0.5 0.5])
        pause(1)
    end
    hold off
end

%% Create array of rays

fprintf('Generating rays source\n')

[P0, k] = CreateRaysArray(illumParam);

if ismember(illumParam.illumCase, [1,3,4])
    P0 = repmat(P0,length(k),1);
    if verbosity == 1 % plotting rays
        hold on,quiver3 (P0(:,3),P0(:,2),P0(:,1),k(:,3)*abs(illumParam.P0(3)),k(:,2)*abs(illumParam.P0(3)),k(:,1)*abs(illumParam.P0(3)),'ShowArrowHead','off','AutoScale','off')
        plot3 (P0(:,3),P0(:,2),P0(:,1),'r*'), hold off
    end
elseif ismember(illumParam.illumCase, [2,5,6])
    k = repmat(k,length(P0),1);
    if verbosity == 1 % plotting rays
        hold on,quiver3 (P0(:,3),P0(:,2),P0(:,1),k(:,3)*abs(illumParam.P0(3)),k(:,2)*abs(illumParam.P0(3)),k(:,1)*abs(illumParam.P0(3)),'ShowArrowHead','off','AutoScale','off')
        plot3 (P0(:,3),P0(:,2),P0(:,1),'r*'), hold off
    end
end

%% Find normals in the surface

fprintf('Surfaces normals calculation\n')

xs = max([abs(min(xu)), max(xu), abs(min(yu)), max(yu)]);                 % x size of grid for interpolation
Nsx = 101;                                                                   % number of divisions in grid for interpolation

xs = linspace(-xs, xs, Nsx);                                           % x axis definition
dxs = xs(2)-xs(1);                                                          % delta size in interpolation grid 

[xs, ys] = meshgrid(xs, xs);                                           % grid creation

zs = cell(1,length(Layer));                                            % interpolated z coordinate per layer/cape 
Normal = cell(1,length(Layer));                                    % cell with and array of normal vector for different positions in the nodes cloud
for l = 1:length(Layer)
    [zs{l}] = griddata(x(Layer{l}),y(Layer{l}),z(Layer{l}),xs,ys);
%     zs{l}(isnan(zs{l})) = 0;
    [Normal{l}] = NormalDir (zs{l}, dxs);
end

% clear zs
%% Ray tracing
thetaref = 0;
[ke1,ke2] = BeamAxes(k', [cosd(thetaref) sind(thetaref) 0]');         % Find polarization axes, based on a linear polarization reference. ke2 is the reference
k11 = -0.5e-12;
k12 = -3.3e-12;
OSC = [k11, k12];
lambda = 532e-9;
n0 = 1.5;                                               % refractive index without load
solMethod = 2;

beamLocNode = zeros(size(k));                % beam location on node
kRefracted = zeros(size(k))';                   % wave vector of refracted beams
beamLoc = cell(1, length(Layer)+1);          % beam location coordinate
beamLoc{1} = P0;
kBeam = cell(1, length(Layer)+1);             % wave vector of light in each layer
kBeam{1} = k';
dn = ones(size(k'));                                % index of refraction in nodes per iteration
dnBeam = cell(1,length(Layer)+2);               % index of refraction in nodes
dnBeam{1} = dn; dnBeam{end} = dn;
ni = ones(length(Layer)+1);                    % index of refraction incident medium
nt = ones(length(Layer)+1);                    % index of refraction transmited medium
beamEllipCoor = cell(length(Layer), length(k));  % Beam direction in ellipsoid coordinate system
nBeam = cell(length(Layer), length(k));             % Index of refraction for a direction of propagation
nBeamDir = cell(length(Layer), length(k));         % Modes of propagation 
globalDCoor = cell(length(Layer), length(k));     % Modes of propagation of D in the global coordinate system
birefringence = cell(1, length(Layer));                % Birefringence
axesRot = cell(1, length(Layer));                         % Axes of rotation of the optical axes
JonesMatrix = cell(1, length(Layer));                   % Jones matrices

for c = 1:length(Layer)                                        % iteration per layer
    fprintf('Ray tracing in layer %u\n',c)
    Strains = cell(1,6);
    for w=1:6                                                        % Stress tensor
    [Strains{w}] = griddata(x(Layer{c}),y(Layer{c}),Data(Layer{c},4+w),xs,ys);
    end
for l = 1: length(k)                                               % iteration per ray
    [beamLocNode(l,:)] =RayPos(x(Layer{c}), y(Layer{c}), z(Layer{c}), beamLoc{c}(l,:), kBeam{c}(:,l));
    [~,CP2B] = sort(sum(abs([xs(:), ys(:)]-beamLocNode(l,1:2)),2)); % CP2B : Closest Point to Beam
    Normal2P = Normal{c}(:,CP2B(1));
    [dn(:,l), StressVD] = StressBir (Strains, beamLocNode(l,1:2), xs(CP2B(1:4)), ys(CP2B(1:4)), CP2B(1:4), OSC, n0);
%     StressVDtmp{l}=StressVD(4:12);
    if c==1
        ni = 1;
        nt = n0+mean(dn(:,l));
    elseif c==length(Layer)
        ni = n0+(mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = 1;
    else
        ni = n0+(mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = n0+mean(dn(:,l));
    end    
    [kRefracted(:,l)] = SnellCalc (kBeam{c}(:,l), Normal2P, ni, nt);   
    [beamEllipCoor{c,l}, Srot, Norder] = SetLocalCord(StressVD(4:6)',StressVD(7:9)',StressVD(10:12)', kRefracted(:,l), ke1(:,l), ke2(:,l));
    [nBeam{c,l}, nBeamDir{c,l}] = IndexEllipCalc(dn(Norder(2),l),dn(Norder(3),l),dn(Norder(1),l),beamEllipCoor{c,l}(:,4), solMethod);
    [globalDCoor{c,l}] = Back2GCS([nBeamDir{c,l}], Srot);
end
clear Strains
beamLoc{c+1} = beamLocNode;
kBeam{c+1} = kRefracted;
dnBeam{c+1} = dn;
[birefringence{c}, axesRot{c}] = JonesMatrixParam(nBeam(c,:), globalDCoor(c,:), kBeam{c+1}, thetaref);
end

fprintf('Jones calculation\n')

for l = 1: length(k)
    Jm = eye(2);
    for c=1:length(Layer)-1
        BirL = (birefringence{c}(l)+birefringence{c+1}(l))/2;
        [Jml] = TwistedBirElem(beamLoc{c+1}(l,:),beamLoc{c+2}(l,:), BirL, lambda, axesRot{c}(l,:),axesRot{c+1}(l,:));
        Jm = Jml*Jm;
    end
    JonesMatrix{l} = Jm;
end

%% Plot polarization map

if verbosity == 1
    In = [1 0]';
    shift = [0 0];
    Efactor = 1;
    OutLayer = length(Layer) ;
    l = 250;
    chiThreshold = pi/125;
    ellipsize = 5;
    arrowsize = 5;
    stepplot = 1:5:51;
    stepplot = sub2ind([51,51], repmat(stepplot,1,length(stepplot)), reshape(repmat(stepplot,length(stepplot),1),1,length(stepplot)^2));
    PosFactor = 5000;
    NumF = 26;
    
    Jones_Ellipse_Plot(JonesMatrix,In,shift, Efactor, NumF, beamLoc{OutLayer}(:,1)*PosFactor, beamLoc{OutLayer}(:,2)*PosFactor, chiThreshold, ellipsize, arrowsize, stepplot)
    hold off
end

% if verbosity == 1
%     In = [1 0]'; In=In/norm(In);
%     shift = [0 0];
%     Efactor = 1;
%     OutLayer = length(Layer) ;
%     l = 250;
%     chiThreshold = pi/125;
%     ellipsize = 5;
%     arrowsize = 5;
% %     stepplot = 1:5:51;
% %     stepplot = sub2ind([51,51], repmat(stepplot,1,length(stepplot)), reshape(repmat(stepplot,length(stepplot),1),1,length(stepplot)^2));
% %     stepplot = 1:3:36;
% %     stepplot = repmat(stepplot,1,15)+36*repmat(0:14,1,length(stepplot));
%     stepplot = 1:1:length(k);
%     PosFactor = 5000;
%     NumF = 28;
%     
%     Jones_Ellipse_Plot(JonesMatrix,In,shift, Efactor, NumF, beamLoc{OutLayer}(:,1)*PosFactor, beamLoc{OutLayer}(:,2)*PosFactor, chiThreshold, ellipsize, arrowsize, stepplot)
%     hold off
% end