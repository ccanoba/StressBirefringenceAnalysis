%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates some optical properties for a stressed optical
% component, employing a ray tracing method. The user must provide
% information regardaing the stresses distribution in the material and the
% illumination comfiguration.
% 
% Inputs :
% 
% Data : Matrix with the information related to the stresses obtained via
%           FEM analysis, the first dimension is the node number. The
%           second dimension are the node index; nodes coordinates x,y,z;
%           normal stresses x,y,z; shear stresses sxy,syz,syz; nodes
%           displacement dx,dy,dz.
% n0 : Refractive index without stress.
% OSC : Optical stress coefficients.
% lambda : Light wavelength.
% illumParam : Illumination parameters defined with SourceDefinition
% Nsx : Number of elements for interpolation.
% thetaref : Reference angle to define the x axis or zero postion of the
%                 model.
% solMethod : Decide whether use the numerical (1) or graphical (2)
%                     approach for the calculation of light eigenmodes and
%                     refractive indices.
% considerDiattenuation : Whether to calculate diattenuation on interfaces
%                                        (1) or not (2).
% verbosity : If (1) plots layers discretization.
% 
% Outputs :
% 
% RetardanceJonesMatrix : Cell array with the retardance Jones matrix for
%                                         every ray in the model.
% Layer : Cell array with the index of every node that belongs to a surface
% beamLoc : Cell array with the coordinates of the ray on each surface of
%                  the model.
% birefringence : Cell array with the two refractive indices percieved by
%                          every ray at the model surface.
% axesRot : Cell array with the orientation of the fast eigenmode of light,
%                 measure from the axis defined with thetaref.
% WF :          Cellarray, the first dimension defines the OPL per ray (1)
%                  or WFE map (2) calculated with the OPL; the second
%                  dimension is the surface, so evolution of the wavefront
%                  can be assessed.
% diattData : Cell array, the first component are the diattenuation matrices
%                  for every ray at the front and rear surfaces; the other
%                  two are the orientation of the perpendicular component
%                  and diattenuation magnitude maps.
% kBeam : Cell array with ray wayvector on every surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [retardanceJonesMatrix, Layer, beamLoc, kBeam, birefringence, axesRot, WF, diattData] = StressBirRayTracing(Data, n0, OSC, lambda, illumParam, Nsx, thetaref, solMethod, considerDiattenuation, verbosity)
% Selection of nodes coordinates from FEM data

xu = Data(:,2); yu = Data(:,3);                                                                     % Undeformed nodes coordinates definition
x = Data(:,2)+Data(:,11); y = Data(:,3)+Data(:,12); z = Data(:,4)+Data(:,13); % Deformed nodes coordinates definition

%% Discretize nodes in the model into layers

fprintf('Nodes cloud discretization\n')

% Variables used to temporarily store the indices of the nodes
% corresponding to each layer.
LayerFound = [];
cont = 1;
haveNodes = true;

% Identify nodes that belongs to a surface in front to light path
while haveNodes == true
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

% Calculate the initial position and wavevector for rays in model
[P0, k] = CreateRaysArray(illumParam);

% Plot rays, Collimated and diverging beams have different initial
% conditions
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
    [zs{l}] = griddata(x(Layer{l}),y(Layer{l}),z(Layer{l}),xs,ys);  % Interpolate surfaces for discrete derivative
    [Normal{l}] = NormalDir (zs{l}, dxs);                       % Calculate normal on every pixel of the interpolated surface
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
retardanceJonesMatrix = cell(1, length(k));       % Jones matrices

for c = 1:length(Layer)                                        % iteration per layer
    fprintf('Ray tracing in layer %u\n',c)
    Strains = cell(1,6);
    for w=1:6                                                        % Stress tensor
    [Strains{w}] = griddata(x(Layer{c}),y(Layer{c}),Data(Layer{c},4+w),xs,ys);  % Interpolates components of the stress tensor
    end
for l = 1: length(k)                                               % iteration per ray
    % Find the ray intersection with the following surface
    [beamLocNode(l,:)] =RayPos(x(Layer{c}), y(Layer{c}), z(Layer{c}), beamLoc{c}(l,:), kBeam{c}(:,l));
    [~,CP2B] = sort(sqrt(sum(([xs(:), ys(:)]-beamLocNode(l,1:2)).^2,2))); % CP2B : Closest Point to Beam
    Normal2P = Normal{c}(:,CP2B(1));
    % Calculates index tensor magnitude and principal directions
    [dn(:,l), StressVD] = StressBir (Strains, beamLocNode(l,1:2), xs(CP2B(1:4)), ys(CP2B(1:4)), CP2B(1:4), OSC, n0);
    % incident and transmitted refractive indices, according to the model surfaces.
    if c==1         % The element is in air
        ni = 1;                                         % index of refraction incident medium
        nt = mean(dn(:,l));                      % index of refraction transmited medium
    elseif c==length(Layer)   % The element is in air
        ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = 1;
    else                                 % Average refractive index
        ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
        nt = mean(dn(:,l));
    end
    % Calculate refracted direction
    [kRefracted(:,l)] = SnellCalc (kBeam{c}(:,l), Normal2P, ni, nt);   
    % Rotate index ellipsoid and wavevector to a local CS defined by index axes
    [beamEllipCoor{c,l}, Srot, Norder] = SetLocalCord(StressVD(4:6)',StressVD(7:9)',StressVD(10:12)', kRefracted(:,l), ke1(:,l), ke2(:,l));
    % Calculates ray eigenvalues and eigenvetors 
    [nBeam{c,l}, nBeamDir{c,l}] = IndexEllipCalc(dn(Norder(2),l),dn(Norder(3),l),dn(Norder(1),l),beamEllipCoor{c,l}(:,4), solMethod);
    % Returns light wavevector and electric field components to global CS
    [globalDCoor{c,l}] = Back2GCS([nBeamDir{c,l}], Srot);
end
clear Strains
% Store relevant parameters in cell
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
        BirL = (birefringence{c}(l)+birefringence{c+1}(l))/2;       % Average retardance between surfaces
        [Jml] = TwistedBirElem(beamLoc{c+1}(l,:),beamLoc{c+2}(l,:), BirL, lambda, axesRot{c}(l,:),axesRot{c+1}(l,:));
        Jm = Jml*Jm;                % Acumulative retardance effect
    end
    retardanceJonesMatrix{l} = Jm;
end
%% Diattenuation calculation at the front and rear faces

diattData = {};
if considerDiattenuation == 1
diattMatrix = cell (length(k),2);   % Diattenuation matrix per ray
diattMag = zeros(length(k),2);    % Diattenuation magnitude
diattAxis = zeros(length(k),2);    % Orientation of transmitance planes

for c=[1, length(Layer)]
    beamLocNode = beamLoc{c+1};
    for l = 1: length(k)
        [~,CP2B] = sort(sum(abs([xs(:), ys(:)]-beamLocNode(l,1:2)),2)); % CP2B : Closest Point to Beam
        Normal2P = Normal{c}(:,CP2B(1));
        if c==1                                 % Calculation on first surface
            ni = 1;
            nt = mean(dn(:,l));
            [diattMatrix{l,1},diattAxis(l,1),diattMag(l,1)] = DiattenuationCalc(kBeam{c}(:,l), kBeam{c+1}(:,l), Normal2P, ni, nt, thetaref);
        elseif c==length(Layer)        % Calculation on last surface
            ni = (mean(dnBeam{c}(:,l))+2*mean(dn(:,l)))/3;
            nt = 1;
            diattMatrix{l,2} = DiattenuationCalc(kBeam{c}(:,l), kBeam{c+1}(:,l), Normal2P, ni, nt, thetaref);
        end
    end
end
diattAxisMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),diattAxis(:,1),xs,ys);
diattMagMap = griddata(beamLoc{2}(:,1),beamLoc{2}(:,2),diattMag(:,1),xs,ys);
diattData = {diattMatrix, diattAxisMap, diattMagMap};
end
%% Wave front error calculation

OPL = zeros(1,length(k));           % optical path length
WF = cell(2,length(Layer)-1);      % Wavefront
for c=1:length(Layer)-1
    for l = 1: length(k)
         if c==1                            % first propagation form source to element does not have birefringence
             n1=0; n2=0;
         else
             % Subtract non stress-induced wavefront error
             n1 = ((nBeam{c-1,l}(1)+nBeam{c,l}(1))/2)-n0;
             n2 = ((nBeam{c-1,l}(2)+nBeam{c,l}(2))/2)-n0;
        end
        OPL(l) = WaveFront(beamLoc{c}(l,:),beamLoc{c+1}(l,:),n1,n2);
    end
    if c==1
        WF{1,c} = OPL;      % Initialize wavefront
    else
        WF{1,c} = WF{1,c-1}+OPL;    % Wavefront error is additive
    end
    WF{2,c} = griddata(beamLoc{c+1}(:,1),beamLoc{c+1}(:,2),WF{1,c},xs,ys);
end