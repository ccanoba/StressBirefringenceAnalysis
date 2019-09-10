%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of nodes defined by its position in space
% and discretizes them into layers. Layers are found assuming the order
% in which a diverging beam would interact with them.
%
% Inputs :
%
% Data : Nodes location and displacement obtained from FEM simulation
% LayerFound : Information of the layer already found
% Param : Illumination source definition
%
% Output :
%
% Layer1 : Node indices for the layer found in this iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Layer1] = FindNodeLayer(Data, LayerFound, Param)

% Displace data that have already been assigned to not take it into
% consideration in further steps.
Data(LayerFound,4) = 1e6;

x = Data(:,2)+Data(:,11);
y = Data(:,3)+Data(:,12);
z = Data(:,4)+Data(:,13);
xu = Data(:,2);
yu = Data(:,3);
zu = Data(:,4);

% Use collimated beam to find the surfaces
if Param.illumCase == 5 || Param.illumCase == 6
    Param.illumCase = Param.illumCase-2;
elseif Param.illumCase == 2
    Param.illumCase = 1;
end

% Define collimated beam
[P0, k] = CreateRaysArray(Param);

%%

distanceSourceNode = sqrt((z-P0(3)).^2);                          % Find distance form source to nodes.
[~,closerNodeSource] = min(distanceSourceNode);          % Find closer node to source.
RP = @(t)P0+[t.*k(:,1),t.*k(:,2),t.*k(:,3)];                               % Position of any point on the ray. RP : Ray position
tmin = (z(closerNodeSource)-P0(3))./k(:,3);                       % Find t value at z for closer node

newCloserNode = RP(tmin);                                              % move point to the closer position, find new closer node

nodeLoc = [];         % nodes location
tR = [];                    % t value in the ray for each nodes

% iterative search of the nodes that first interact with the rays coming
% from the source

for m = 1:length(newCloserNode)
    searchNodeFlag = true;               % Search for node. SN : Search Node.
    nodeLocFlag = closerNodeSource;       % closer node in previous iteration. ndrefp : closer node position.
    
    while searchNodeFlag == true
        distanceSourceNode = sqrt((x-newCloserNode(m,1)).^2+(y-newCloserNode(m,2)).^2+(z-newCloserNode(m,3)).^2);
        [~,currentNodeLoc] = min(distanceSourceNode);
        tmin = (z(currentNodeLoc)-P0(3))./k(m,3);
        if currentNodeLoc == nodeLocFlag
            searchNodeFlag = false;
            nodeLoc = [nodeLoc currentNodeLoc];
            tR = [tR tmin];
        end
        nodeLocFlag = currentNodeLoc;
    end
end

% fit nodes positions into a quadratic surface. Generated with curvefitting-matlab
[Scoeff, checktmp] = FitSurface(xu(nodeLoc),yu(nodeLoc),zu(nodeLoc));

zL = @(x,y) Scoeff.p00 + Scoeff.p10.*x + Scoeff.p01.*y + Scoeff.p20.*x.^2 + Scoeff.p11.*x.*y + Scoeff.p02.*y.^2;

% Find z position for each node depending on its x-y coordinates
ZL = zL(xu,yu);               % With xu,yu,zu

% Threshold to determine which nodes belongs to the surface/ Sensitive step
if checktmp.rmse>1e-10
[Layer1] = find(abs(zu-ZL)<=4.5*checktmp.rmse);        % With xu,yu,zu
else
[Layer1] = find(abs(zu-ZL)<=10*max(abs(zu(nodeLoc)-ZL(nodeLoc))));        % With xu,yu,zu
end