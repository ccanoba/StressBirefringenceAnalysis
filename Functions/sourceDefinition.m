%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function generates the struct file with the parameters required for
%  the light source definition.

function [sourceParam] = sourceDefinition(varargin)

sourceParam.RI = [];
sourceParam.NRI = [];
sourceParam.NTI = [];
sourceParam.P0 = [];
sourceParam.posLim = [];
sourceParam.xside = [];
sourceParam.yside = [];
sourceParam.Nx = [];
sourceParam.Ny = [];
sourceParam.illumCase = [];

if varargin{1} == 1 || varargin{1} == 2 % Diverging spherical discretization or Collimated spherical discretization
    sourceParam.illumCase = varargin{1};
    sourceParam.RI = varargin{2};
    sourceParam.NRI = varargin{3};
    sourceParam.NTI = varargin{4};
    sourceParam.P0 = varargin{5};
    sourceParam.posLim = varargin{6};
    
elseif varargin{1} == 3 || varargin{1} == 5 % Diverging square discretization. Square shape or Collimated square discretization. Square shape
    
    sourceParam.illumCase = varargin{1};
    sourceParam.xside = varargin{2};
    sourceParam.yside = varargin{3};
    sourceParam.Nx = varargin{4};
    sourceParam.Ny = varargin{5};
    sourceParam.P0 = varargin{6};
    sourceParam.posLim = varargin{7};
    
elseif varargin{1} == 4 || varargin{1} == 6 % Diverging square discretization. Circular shape or Collimated square discretization. Circular shape
    
    sourceParam.illumCase = varargin{1};
    sourceParam.RI = varargin{2};
    sourceParam.xside = varargin{3};
    sourceParam.yside = varargin{4};
    sourceParam.Nx = varargin{5};
    sourceParam.Ny = varargin{6};
    sourceParam.P0 = varargin{7};
    sourceParam.posLim = varargin{8};
    
else
    
    fprintf('That case does not exist\n')
    
end

