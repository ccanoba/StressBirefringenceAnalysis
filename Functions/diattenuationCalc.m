function [diattMatrix, gammad, diattenuation] = diattenuationCalc(Kin, Kref, N, n1, n2, thetaref)

ts = (2*n1*dot(N,Kin))/(n1*dot(N,Kin)+n2*dot(N,Kref));
tp = (2*n1*dot(N,Kin))/(n2*dot(N,Kin)+n1*dot(N,Kref));

[~,ke2] = BeamAxes(Kin, [cosd(thetaref) sind(thetaref) 0]'); 

Vs = cross(Kin, N);
Vs = Vs/norm(Vs);
% Vt = cross(Kin, Vs);

ST = cross(Vs,ke2);
STsign = sign(dot(ST,Kin));
ST = sqrt(ST(1).^2+ST(2).^2+ST(3).^2);
CT = dot(Vs,ke2);
CT = CT.*STsign;

gammad = atan2d(ST,CT)';
gammad = (gammad-90);

% gammad = asin(norm(cross(Vs, ke2)));
if isnan(gammad)
    gammad = 0;
end

diattenuation = abs(ts.^2-tp.^2)/(ts.^2+tp.^2);

diattMatrix = [ts*(cosd(gammad))^2+tp*(sind(gammad))^2 (ts-tp)*cosd(gammad)*sind(gammad);...
                      (ts-tp)*cosd(gammad)*sind(gammad) ts*(sind(gammad))^2+tp*(cosd(gammad))^2];