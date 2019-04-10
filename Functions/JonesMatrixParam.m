%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This function calculates the birefringence and fast axis angle respect
%  to the beam coordinate system
% 
% inputs
% 
% Nbeam : index of refraction of the medium
% GlobalDCoor : modes of vibration of light in the medium
% kL : direction of propagation of light into the medium
% 
% outputs
% 
% birefringence : birefringence in the evaluated layer
% axesRot : angle between mode of vibration and beam reference plane
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [birefringence, axesRot] = JonesMatrixParam(Nbeam, GlobalDCoor, kL, thetaref)

dn = Nbeam;
dn = real(cell2mat(dn));
dn1=dn(1:2:end);
dn2=dn(2:2:end);
birefringence = dn1-dn2;

D=GlobalDCoor;
D=real(cell2mat(D));
D2=D(:,2:2:end);

[~,ke2] = BeamAxes(kL, [cosd(thetaref) sind(thetaref) 0]'); 
ST = cross(D2,ke2);
STsign = sign(dot(ST,kL));
ST = sqrt(ST(1,:).^2+ST(2,:).^2+ST(3,:).^2);
CT = dot(D2,ke2);
CT = CT.*STsign;

axesRot = atan2d(ST,CT)';
axesRot = axesRot-90;
end