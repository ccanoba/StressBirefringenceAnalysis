function [psi,chi]=Jones2angle(Jones)

Ax=sqrt(Jones(1,1)*conj(Jones(1,1)));
Ay=sqrt(Jones(2,1)*conj(Jones(2,1)));
dx=atan2(imag(Jones(1,1)),real(Jones(1,1)));
dy=atan2(imag(Jones(2,1)),real(Jones(2,1)));
d=dy-dx;
psi=1/2*atan2((2*Ax*Ay)*cos(d),(Ax^2-Ay^2));
chi=1/2*asin(((2*Ax*Ay)/(Ax^2+Ay^2))*sin(d));