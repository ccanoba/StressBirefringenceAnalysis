%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes a set of nodes defined by its position in space
% and discretizes them into layers. Layers are found assuming the order
% in which a diverging beam would interact with them.
%
% 


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

minz=min(z);

% parameter for the definition of illumination

% RI = Param.RI;
% % Number of divisions in RI
% NRI = Param.NRI;
% % Number of divisions in theta
% NTI = Param.NTI;
% % Distance to illumination plane
% ZI = abs(Param.P0(3))+Param.posLim(5);
% % Source Point
% P0 = Param.P0;
% 
% % RI = 12e-3;                % Illuminaion cone radio
% % NRI = 15;                   % Number of divisions in RI
% % NTI = 36;                   % Number of divisions in theta
% % ZI = 200e-3;              % Distance to illumination plane
% % ZI = ZI+minz;
% % 
% % P0 = [0 0e-3 -200e-3]; % Source Point Location
% 
% % Create rays wave vectors
% % VRI : radial vector component
% % VTI : angular vector component
% VRI = RI/NRI:RI/NRI:RI;
% VRI = repmat(VRI,NTI,1);
% VRI=reshape(VRI,[1,numel(VRI)]);
% VTI = 0: 2*pi/NTI: 2*pi-(2*pi/NTI);
% VTI=repmat(VTI,1,NRI);
% % k : wavevector
% k = [VRI'.*cos(VTI'), VRI'.*sin(VTI'), ones(NRI*NTI,1)*ZI];
% k = [[0 0 ZI]; k];
% k = k./sqrt(k(:,1).^2+k(:,2).^2+k(:,3).^2);  % normalization of k

if Param.illumCase == 5 || Param.illumCase == 6
    Param.illumCase = Param.illumCase-2;
elseif Param.illumCase == 2
    Param.illumCase = 1;
end

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

% Threshold to determine which nodes belongs to the surface
if checktmp.rmse>1e-10
[Layer1] = find(abs(zu-ZL)<=4.5*checktmp.rmse);        % With xu,yu,zu
else
[Layer1] = find(abs(zu-ZL)<=10*max(abs(zu(nodeLoc)-ZL(nodeLoc))));        % With xu,yu,zu
end