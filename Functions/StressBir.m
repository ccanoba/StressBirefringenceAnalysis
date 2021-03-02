%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This function calculates the changes in the index of refraction of a
% material due to photoelactic effect.
%
% Inputs
% 
% Strains : interpolated components of the stress tensor
% Pos : beam coordinates
% x-y : x-y grid coordinates
% CN2P : index of the closest node to beam
% OSC : Optical Stress Coefficient
% n0 : index of refraction of the material without load
% 
% Output
% 
% dn : changes in the index of refraction due to stress birefringence
% StressVD : principal directions of the index ellipsoid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors:  Camilo Cano {1*}, Pablo Zuluaga-Ram�rez {2}, Ren� Restrepo {1,3}
%   1. Applied Optics Group, Universidad EAFIT, Carrera 49 # 7 Sur-50,
%   Medell�n, Colombia.
%   2. European Southern Observatory Headquarters, Karl-Schwarzschild-Str. 2, 
%   85748 Garching bei Munchen, Germany
%   3. Aerospace Optics Instrumentation Division, National Institute of Aerospace
%   Technology - INTA, Ctra de Ajalvir, Km 4, Torrejon de Ardoz, 28850 Madrid, Spain
%	* <ccanoba@eafit.edu.co>    -   2019

function [dn, StressVD] = StressBir (Strains, Pos, x, y, CN2P, OSC, n0)

% Calculates index ellipsoid due to stresses

StressT = zeros(1,6);

% bilinear interpolation to find stresses where the ray fall.

for i=1:6
    StressT(i) = InterPolStress(Strains{i}, Pos, x, y, CN2P);
end

stressTensor = [StressT(1), StressT(4), StressT(6);...
                          StressT(4), StressT(2), StressT(5);...
                          StressT(6), StressT(5), StressT(3)];
                      
% Diagonalization of the stress tensor                      
[PSD, PSV]=principalStresses(stressTensor);

StressVD(1:3) = PSV([1,5,9]);
StressVD(4:12) = PSD(1:9); 

% made all stress direntions point in the positive z axis.
if StressVD(6)<0
    StressVD(4:6)=StressVD(4:6)*-1;
end
if StressVD(9)<0
    StressVD(7:9)=StressVD(7:9)*-1;
end
if StressVD(12)<0
    StressVD(10:12)=StressVD(10:12)*-1;
end

dn1 = OSC(1)*StressVD(:,1)+OSC(2)*(StressVD(:,2)+StressVD(:,3));
dn2 = OSC(1)*StressVD(:,2)+OSC(2)*(StressVD(:,3)+StressVD(:,1));
dn3 = OSC(1)*StressVD(:,3)+OSC(2)*(StressVD(:,1)+StressVD(:,2));

dn = n0+[dn1, dn2, dn3];

    function [Fv] = InterPolStress(Strains, Pos, x, y, CN2P)
%       Bilinear interpolation function
% 
%       inputs
%       Strains : Data to interpolate
%       Pos : coordinates in which data will be interpolated
%       x-y : x-y coordinates of data
%       CN2P : index of the data close to the interpolation region
%
%       output
% 
%       Fv : interpolated value

        Strains(isnan(Strains)) = 0;

        xs = sort(x);
        ys = sort(y);
        
        if sum(x==Pos(1)&y==Pos(2))~=0
            Fv = Strains(CN2P(x==Pos(1)&y==Pos(2)));
        elseif (xs(1)~=xs(2)||xs(3)~=xs(4)||ys(1)~=ys(2)||ys(3)~=ys(4))
            [~,vectcomp]=max([length(unique(x)),length(unique(y))]);
            index = 1:4;
            if vectcomp ==1
                [~,uniquex] = unique(x);
                RepElem = index(~ismember(index,uniquex));
                RepElem = ismember(x,x(RepElem));
                y = y(RepElem);
                Fv = (max(y)-Pos(2))/(max(y)-min(y))*Strains(CN2P(y==min(y)))+(Pos(2)-min(y))/(max(y)-min(y))*Strains(CN2P(y==max(y)));
            else
                [~,uniquey] = unique(y);
                RepElem = index(~ismember(index,uniquey));
                RepElem = ismember(y,y(RepElem));
                x = x(RepElem);
                Fv = (max(x)-Pos(1))/(max(x)-min(x))*Strains(CN2P(x==min(x)))+(Pos(1)-min(x))/(max(x)-min(x))*Strains(CN2P(x==max(x)));
            end
        else
            
            x1 = min(x); x2 = max(x);
            y1 = min(y); y2 = max(y);
            
            Q11=Strains(CN2P(x==x1 & y==y1));
            Q12=Strains(CN2P(x==x1 & y==y2));
            Q21=Strains(CN2P(x==x2 & y==y1));
            Q22=Strains(CN2P(x==x2 & y==y2));
            
            Fv = (1/((x2-x1)*(y2-y1)))*[x2-Pos(1) Pos(1)-x1]*[Q11 Q12; Q21 Q22]*[y2-Pos(2);Pos(2)-y1];
        end
    end
end