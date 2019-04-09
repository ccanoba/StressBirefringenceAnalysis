function [axes1, axes2] = BeamAxes(k, Refaxis)

% axes of the beam components, considering and input plane given by the
% Refaxis,

% Refaxis = [1 0 0]';
axes1 = cross(k,repmat(Refaxis,1,size(k,2)));
axes1 = axes1./sqrt(axes1(1,:).^2+axes1(2,:).^2+axes1(3,:).^2); 
axes2 = cross(k,axes1);
axes2 = axes2./sqrt(axes2(1,:).^2+axes2(2,:).^2+axes2(3,:).^2); 

end



