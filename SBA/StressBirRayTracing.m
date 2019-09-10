%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs :
% 
% Data:
% n0 :
% OSC : 
% lambda :
% illumParam :
% Nsx :
% thetaref :
% solMethod :
% considerDiattenuation:
% verbosity :

% Outputs :
% 
% RetardanceJonesMatrix :
% Layer :
% beamLoc :
% birefringence :
% axesRot :
% WF :
% diattData :
% kBeam :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RetardanceJonesMatrix, Layer, beamLoc, birefringence, axesRot, WF, diattData, kBeam] = StressBirRayTracing(Data, n0, OSC, lambda, illumParam, Nsx, thetaref, solMethod, considerDiattenuation, verbosity)

xu = Data(:,2); yu = Data(:,3);                                                                     % Undeformed nodes coordinates definition
x = Data(:,2)+Data(:,11); y = Data(:,3)+Data(:,12); z = Data(:,4)+Data(:,13); % Deformed nodes coordinates definition

%% Discretize nodes in the model into layers

fprintf('Nodes cloud discretization\n')

% Variables used to temporarily store the indices of the nodes
% corresponding to each layer.
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

xs = max([abs(min(xu)), max(xu), abs(min(yu)), max(yu)]);                 % size of grid for interpolation

xs = linspace(-xs, xs, Nsx);                                           % axis definition
dxs = xs(2)-xs(1);                                                          % delta size in interpolation grid 

[xs, ys] = meshgrid(xs, xs);                                           % grid creation

zs = cell(1,length(Layer));                                            % interpolated z coordinate per layer/cape 
Normal = cell(1,length(Layer));                                    % cell with and array of normal vector for different positions in the nodes cloud
for l = 1:length(Layer)
    [zs{l}] = griddata(x(Layer{l}),y(Layer{l}),z(Layer{l}),xs,ys);
    [Normal{l}] = NormalDir (zs{l}, dxs);
end
clear zs
%% Ray tracing

[ke1,ke2] = BeamAxes(k', [cosd(thetaref) sind(thetaref) 0]');         % Find polarization axes, based on a linear polarization reference. ke2 is the reference
beamLocNode = zeros(size(k));                % beam location on node
kRefracted = zeros(size(k))';                   % wave vector of refracted beams
beamLoc = cell(1, length(Layer)+1);          % beam location coordinate
beamLoc{1} = P0;
kBeam = cell(1, length(Layer)+1);             % wave vector of light in each layer
kBeam{1} = k';
dn = ones(size(k'));                                % index of refraction in nodes per iteration
dnBeam = cell(1,length(Layer)+2);               % index of refraction in nodes
dnBeam{1} = dn; dnBeam{end} = dn;
beamEllipCoor = cell(length(Layer), length(k));  % Beam direction in ellipsoid coordinate system
nBeam = cell(length(Layer), length(k));             % Index of refraction for a direction of propagation
nBeamDir = cell(length(Layer), length(k));         % Modes of propagation 
globalDCoor = cell(length(Layer), length(k));     % Modes of propagation of D in the global coordinate system
birefringence = cell(1, length(Layer));                % Birefringence
axesRot = cell(1, length(Layer));                         % Axes of rotation of the optical axes
RetardanceJonesMatrix = cell(1, length(k));                   % Jones matrices

for c = 1:length(Layer)                                        % iteration per layer
    fprintf('Ray tracing in layer %u\n',c)
    Strains = cell(1,6);
    for w=1:6                                                        % Stress tensor
    [Strains{w}] = griddata(x(Layer{c}),y(Layer{c}),Data(Layer{c},4+w),xs,ys);
    end
for l = 1: length(k)                                               % iteration per ray
    [beamLocNode(l,:)] =RayPos(x(Layer{c}), y(Layer{c}), z(Layer{c}), beamLoc{c}(l,:), kBeam{c}(:,l));
    [~,CP2B] = sort(sqrt(sum(([xs(:), ys(:)]-beamLocNode(l,1:2)).^2,2))); % CP2B : Closest Point to Beam
    Normal2P = Normal{c}(:,CP2B(1));
    [dn(:,l), StressVD] = StressBir (Strains, beamLocNode(l,1:2), xs(CP2B(1:4)), ys(CP2B(1:4)), CP2B(1:4), OSC, n0);
    if c==1
        ni = 1;                                         % index of refraction incident medium
        nt = mean(dn(:,l));                      % index of refraction transmited medium
    elseif c==length(Layer)
        ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = 1;
    else
        ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = mean(dn(:,l));
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

%% Calculation of retardance Jones matrices, assuming twisted nematic behaviour of the medium
fprintf('Jones calculation\n')

for l = 1: length(k)
    Jm = eye(2);
    for c=1:length(Layer)-1
        BirL = (birefringence{c}(l)+birefringence{c+1}(l))/2;
        [Jml] = TwistedBirElem(beamLoc{c+1}(l,:),beamLoc{c+2}(l,:), BirL, lambda, axesRot{c}(l,:),axesRot{c+1}(l,:));
        Jm = Jml*Jm;
    end
    RetardanceJonesMatrix{l} = Jm;
end
%% Diattenuation calculation at the front and rear faces

diattData = {};
if considerDiattenuation == 1
diattMatrix = cell (length(k),2);
diattMag = zeros(length(k),2);
diattAxis = zeros(length(k),2);

for c=[1, length(Layer)]
    beamLocNode = beamLoc{c+1};
    for l = 1: length(k)
        [~,CP2B] = sort(sum(abs([xs(:), ys(:)]-beamLocNode(l,1:2)),2)); % CP2B : Closest Point to Beam
        Normal2P = Normal{c}(:,CP2B(1));
        if c==1
            ni = 1;
            nt = mean(dn(:,l));
            [diattMatrix{l,1},diattAxis(l,1),diattMag(l,1)] = diattenuationCalc(kBeam{c}(:,l), kBeam{c+1}(:,l), Normal2P, ni, nt, thetaref);
        elseif c==length(Layer)
            ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
            nt = 1;
            diattMatrix{l,2} = diattenuationCalc(kBeam{c}(:,l), kBeam{c+1}(:,l), Normal2P, ni, nt, thetaref);
        end
    end
end
diattAxisMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),diattAxis(:,1),xs,ys);
diattMagMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),diattMag(:,1),xs,ys);
diattData = {diattMatrix, diattAxisMap, diattMagMap};
end
%% Wave front error calculation

OPL = zeros(1,length(k));
WF = cell(2,length(Layer)-1); 
for c=1:length(Layer)-1
    for l = 1: length(k)
         if c==1
             n1=0; n2=0;
         else
             n1 = ((nBeam{c-1,l}(1)+nBeam{c,l}(1))/2)-n0;
             n2 = ((nBeam{c-1,l}(2)+nBeam{c,l}(2))/2)-n0;
        end
        OPL(l) = waveFront(beamLoc{c}(l,:),beamLoc{c+1}(l,:),n1,n2);
    end
    if c==1
        WF{1,c} = OPL;
    else
        WF{1,c} = WF{1,c-1}+OPL;
    end
    WF{2,c} = griddata(beamLoc{c+1}(:,1),beamLoc{c+1}(:,2),WF{1,c},xs,ys);
end