%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function determines the position in which will lay a ray that
% propagates from a origin point P0 in a direction k, over a surface
% defined by a set of nodes, assuming that those represent a
% discrete representation of a surface.
%
% Inputs
%
% x, y, z : cartesian coordinates of the nodes that describe the surface
% P0 : Origin of the beam. (x,y,z coordinates)
% k : wave vector, direction of propagation of the beam. (vector in cartesian coordinates)
%
% Output
%
% beamPos : Cartesian coordinates of the beam over the surface
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

function [beamPos] = RayPos(x, y, z, P0, k)

[d1, CN1, dNB] = Point2LineDist(x, y, z, P0, k);

dNB(CN1) = inf;  

for i = 1:8         % Find 8 other closest nodes
    
    [d2(i),CN2(i)]=min(dNB);
    
    PC1 = [x(CN1), y(CN1), z(CN1)];
    PC2 = [x(CN2(i)), y(CN2(i)), z(CN2(i))];
    kL(i,:) = PC2-PC1;
    
    Tc = (((d1-d2(i))/(d1+d2(i)))+1)/2;
    PObj(i,:) = PC1+kL(i,:)*Tc;
    dNB(CN2(i)) = inf;
end

[~, CNObj] = Point2LineDist(PObj(:,1), PObj(:,2), PObj(:,3), P0, k);

beamPos = PObj(CNObj,:);

    function [d, CN, dPr] = Point2LineDist(x, y, z, P0, k)
        
%         This function determines the closest point on a line respect to a
%         cloud of points
%         Equation is described in https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
%         Inputs
%
%         x, y, z : cartesian coordinates of the nodes that describe the surface
%         P0 : Origin of the beam. (x,y,z coordinates)
%         k : wave vector, direction of propagation of the beam. (vector in cartesian coordinates)
%         Outputs
%
%         dPr : vector with de closest distance from every node to the ray
%         CN : index of the closest node to the ray
%         d : distance from the closest node to the ray
%
%         Variables notation are acording to matemathical description
        
        AP = P0-[x, y, z];
        
        kxAP = cross(repmat(k,1,length(AP)), AP');
        
        NkxAP = sqrt(kxAP(1,:).^2+kxAP(2,:).^2+kxAP(3,:).^2);
        
        Nk = norm(k);
        
        dPr = NkxAP/Nk;
        
        [d,CN]=min(dPr);
    end
end