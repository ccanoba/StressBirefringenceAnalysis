function [OPL] = waveFront(Pi, Pf, n1, n2) 

L = sqrt((Pf(:,1)-Pi(:,1)).^2+(Pf(:,2)-Pi(:,2)).^2+(Pf(:,3)-Pi(:,3)).^2);   % Lenght of the medium
OPL = L*(n1+n2)/(2);
