%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates a rays array that will interact with the volume
% described by the nodes of the finite element model
%
% Input :
%
% sourceParam : The source parameter are defined using sourceDefinition
%
% Output :
%
% P0 : Origin of every ray in the model
% k : wavevector for every ray in the model
%
% case 1 : Diverging spherical discretization
% case 2 : Collimated spherical discretization
% case 3 : Diverging square discretization. Square shape
% case 4 : Diverging square discretization. Circular shape
% case 5 : Collimated square discretization. Square shape
% case 6 : Collimated square discretization. Circular shape
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

function [P0, k] = CreateRaysArray(sourceParam)

illumCase = sourceParam.illumCase;

if illumCase == 1 || illumCase==2                % Spherical segmentation. Theta and R
    % Illuminaion cone radio
    RI = sourceParam.RI;
    % Number of divisions in RI
    NRI = sourceParam.NRI;
    % Number of divisions in theta
    NTI = sourceParam.NTI;
    % Distance to illumination plane
    ZI = abs(sourceParam.P0(3))+sourceParam.posLim(5);
    % Source Point
    P0 = sourceParam.P0;
    
    % Create rays wave vectors
    VRI = RI/NRI:RI/NRI:RI;
    VRI = repmat(VRI,NTI,1);
    VRI=reshape(VRI,[1,numel(VRI)]);
    VTI = 0: 2*pi/NTI: 2*pi-(2*pi/NTI);
    VTI=repmat(VTI,1,NRI);    
    % Definition of a point diveging source with circular shape and polar
    % discretization
    if illumCase == 1
        k = [VRI'.*cos(VTI'), VRI'.*sin(VTI'), ones(NRI*NTI,1)*ZI];
        k = [[0 0 ZI]; k];
        k = k./sqrt(k(:,1).^2+k(:,2).^2+k(:,3).^2);  % normalization of k
    else
        % Definition of a collimated illumination with circular shape and
        % polar discretization
        P0 = [VRI'.*cos(VTI'), VRI'.*sin(VTI'), sign(sourceParam.P0(3))*ones(numel(VRI),1)*ZI];
        k = [0,0,1];
    end    
    
elseif illumCase == 3 || illumCase == 4
    % Square discretization, square aperture
    xside = sourceParam.xside;
    yside = sourceParam.yside;
    Nx = sourceParam.Nx;
    Ny = sourceParam.Ny;
   
    ZI = abs(sourceParam.P0(3))+sourceParam.posLim(5);
    P0 = sourceParam.P0; % Source Point
    
    x = linspace(-xside/2, xside/2, Nx+1);
    y = linspace(-yside/2, yside/2, Ny+1);
    
    [x,y] = meshgrid(x,y);  % Square discretization, square aperture
    
    if illumCase == 4 % Square discretization, circular aperture
        [~,R] = cart2pol(x,y);
        
        mask = find(R>sourceParam.RI);
        
        x(mask) = [];
        y(mask) = [];
    end
    
    k = [x(:),y(:), ones(numel(x),1)*ZI];
    
    k = k./sqrt(k(:,1).^2+k(:,2).^2+k(:,3).^2);  % normalization of k
    
elseif illumCase == 5 || illumCase == 6 % Collimated field
    
    xside = sourceParam.xside;
    yside = sourceParam.yside;
    Nx = sourceParam.Nx;
    Ny = sourceParam.Ny;
    ZI = abs(sourceParam.P0(3))+sourceParam.posLim(5);
    
    x = linspace(-xside/2, xside/2, Nx+1);
    y = linspace(-yside/2, yside/2, Ny+1);
    
    [x,y] = meshgrid(x,y);
    
    if illumCase == 5 % Square aperture
        
        P0 = [x(:), y(:), sign(sourceParam.P0(3))*ones(numel(x),1)*ZI];
        
    elseif illumCase == 6 % Circular aperture
        
        [~,R] = cart2pol(x,y);
        
        mask = find(R<sourceParam.RI);
        
        P0 = [x(mask), y(mask), sign(sourceParam.P0(3))*ones(length(mask),1)*ZI];
        
    end
    
    k = [0,0, ZI];
    
    k = k./sqrt(k(:,1).^2+k(:,2).^2+k(:,3).^2);  % normalization of k
    
else
    error (sprintf('Error. \nThat case does not exist'))
end
end
