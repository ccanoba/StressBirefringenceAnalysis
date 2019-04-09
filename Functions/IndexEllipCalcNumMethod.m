%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Numerical method to calculate efective birefringence and modes of
% vibration for a beam propagating in an anisotropic media.
% 
% inputs
% 
% nx,ny,nz : components of the index ellipsoid
% S : wavevector of light
%
% outputs
% 
% Nbeam : indices of refraction for the components of E
% NbeamDir : directions of vibration of E
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Nbeam, NbeamDir] = IndexEllipCalcNumMethod(nx,ny,nz,S)

ux = S(1).^2;
uy = S(2).^2;
uz = S(3).^2;
ex = (nx)^2;
ey = (ny)^2;
ez = (nz)^2;

A = ex*ux+ey*uy+ez*uz;
B = ux*ex*(ey+ez)+uy*ey*(ex+ez)+uz*ez*(ex+ey);
C = ex*ey*ez;

N1E = B/(2*A)+sqrt(B^2-4*A*C)/(2*A);
N2E = B/(2*A)-sqrt(B^2-4*A*C)/(2*A);

Nbeam = [sqrt(N1E), sqrt(N2E)];
[Nbeam, orderN] = sort(Nbeam);

Ev1 = [S(1)/(N1E-ex), S(2)/(N1E-ey), S(3)/(N1E-ez)]';
Ev2 = [S(1)/(N2E-ex), S(2)/(N2E-ey), S(3)/(N2E-ez)]';

if mean(isnan(Ev1))~=0 || mean(isnan(Ev2))~=0
    Ev1(isnan(Ev1))=1e10;
    Ev2(isnan(Ev2))=1e10;
end

Ev1 = Ev1/norm(Ev1);
Ev2 = Ev2/norm(Ev2);

if abs(Ev1'*S)>0.5 
    Ev1 = cross(Ev2,S);
    Ev1 = Ev1/norm(Ev1);
elseif abs(Ev2'*S)>0.5
    Ev2 = cross(Ev1,S);
    Ev2 = Ev2/norm(Ev2);
end


NbeamDir = [Ev1, Ev2];
NbeamDir(:, orderN) = NbeamDir;

[~,maxvx] = max(abs(NbeamDir(1,:)));
[~,maxvy] = max(abs(NbeamDir(2,:)));

if NbeamDir(1, maxvx)<0
    NbeamDir(:,maxvx)=NbeamDir(:,maxvx)*-1;
end
if NbeamDir(2, maxvy)<0
    NbeamDir(:,maxvy)=NbeamDir(:,maxvy)*-1;
end
end