%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Graphical method to calculate efective birefringence and modes of
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
% OpAxis1 : optical axis direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Nbeam, NbeamDir, OpAxis1] = IndexEllipCalcGraphMethod(nx,ny,nz,S)

% To find optic axes
ka = @(na,nb,nc) (nc.*sqrt((nb.^2-na.^2)./(nc.^2-na.^2)));
kb = @(nb,ka) sqrt(nb+ka).*sqrt(nb-ka);
Ellipeq = @(v) v(1)^2/nx^2+v(2)^2/ny^2+v(3)^2/nz^2;

% Propose possible permutations
na = [nz nx ny];
nb = [nx ny nz];
nc = [ny nz nx];

% Calculate axes
pa = ka(na,nb,nc);
pb = kb(nb,pa);

% Find axes orientation
RotAxes = find(real(pa)~=0 & real(pb)~=0);
if isempty(RotAxes) || length(RotAxes)>1
    RotAxes=4;
end
switch RotAxes
    case 1
        OpAxis1 = [0 pb(RotAxes) pa(RotAxes)];
        OpAxis2 = [0 pb(RotAxes) -pa(RotAxes)];
        OpAxis1 = OpAxis1/norm(OpAxis1);
        OpAxis2 = OpAxis2/norm(OpAxis2);
    case 2
        OpAxis1 = [pa(RotAxes) 0 pb(RotAxes)];
        OpAxis2 = [-pa(RotAxes) 0 pb(RotAxes)];
        OpAxis1 = OpAxis1/norm(OpAxis1);
        OpAxis2 = OpAxis2/norm(OpAxis2);
    case 3
        OpAxis1 = [pb(RotAxes) pa(RotAxes) 0];
        OpAxis2 = [pb(RotAxes) -pa(RotAxes) 0];
        OpAxis1 = OpAxis1/norm(OpAxis1);
        OpAxis2 = OpAxis2/norm(OpAxis2);
    otherwise
        Nbeam = [nx, ny];
        NbeamDir = [[1 0 0]', [0 1 0]'];
        OpAxis1 = [0 0 0];
        return
end

[Aye,Bye]=EllipsoidPlaneIntersection(S(1),S(2),S(3),0,nx,ny,nz);

R1 = cross(S, OpAxis1);
R1 = R1/norm(R1);
R2 = cross(S, OpAxis2);
R2 = R2/norm(R2);

R3 = R1+R2;
if norm(R3) == 0
    R2 = -1*R2;
    R3 = R1+R2;
end
R3 = R3/norm(R3);

R4 = cross(R3,S);
R4 = R4/norm(R4);

Nbeam = [Aye,Bye];
[~,NDirCor]=min(abs([Ellipeq(R3*Aye), Ellipeq(R3*Bye)]-1));
if NDirCor == 1
    NbeamDir = [R3', R4'];
else
    NbeamDir = [R4', R3'];
end

[~,Nmag]=sort(Nbeam);

Nbeam = Nbeam(Nmag);
NbeamDir = NbeamDir(:,Nmag);

[~,maxvx] = max(abs(NbeamDir(1,:)));
[~,maxvy] = max(abs(NbeamDir(2,:)));

if NbeamDir(1, maxvx)<0
    NbeamDir(:,maxvx)=NbeamDir(:,maxvx)*-1;
end
if NbeamDir(2, maxvy)<0
    NbeamDir(:,maxvy)=NbeamDir(:,maxvy)*-1;
end