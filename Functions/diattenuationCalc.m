function [diattMatrix, gammad, diattenuation] = diattenuationCalc(Kin, Kref, N, n1, n2, thetaref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This funciton calculates the diatenuation matrix on an interface based on
% Fresnel's equations
%
% Inputs :
%
% Kin : Incident direction of propagation defined by a unit vector
% Kref : Refracted direction of propagation defined by a unit vector
% N : Normal surface direction  defined by a unit vector
% n1 : Incidence refractive index
% n2 : Refracted refractive index
% thetaref : Angle used to define polarized plane, 0 is the x-axis of
%                global CS
%
% Output :
%
% diattMatrix : Diattenuation matrix for every ray in the model
% gammad : Matrix with the orientation angle for the eigenmodes of light
%                   over the system aperture
% diattenuation : Matrix with the diattenuation value over the system aperture
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perpendicular and parallel transmitances calculation, using Fresnel's
% equations
ts = (2*n1*dot(N,Kin))/(n1*dot(N,Kin)+n2*dot(N,Kref));
tp = (2*n1*dot(N,Kin))/(n2*dot(N,Kin)+n1*dot(N,Kref));

% Calculation of the orientation diattenuation planes
[~,ke2] = BeamAxes(Kin, [cosd(thetaref) sind(thetaref) 0]'); 
Vs = cross(Kin, N);
Vs = Vs/norm(Vs);
ST = cross(Vs,ke2);
STsign = sign(dot(ST,Kin));
ST = sqrt(ST(1).^2+ST(2).^2+ST(3).^2);
CT = dot(Vs,ke2);
CT = CT.*STsign;
gammad = atan2d(ST,CT)';
gammad = (gammad-90);
if isnan(gammad)
    gammad = 0;
end

% Diattenuation calculation
diattenuation = abs(ts.^2-tp.^2)/(ts.^2+tp.^2);

% Diattenuation matrix definition
diattMatrix = [ts*(cosd(gammad))^2+tp*(sind(gammad))^2 (ts-tp)*cosd(gammad)*sind(gammad);...
                      (ts-tp)*cosd(gammad)*sind(gammad) ts*(sind(gammad))^2+tp*(cosd(gammad))^2];