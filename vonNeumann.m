function [stability] = vonNeumann(dt, dx, alpha)
%% run von Neumann stability analysis for 1D heat equation
% inputs:
%   dt = timestep (s)
%   dx = spatial step (m)
%   alpha = vector of thermal diffusivities (m2/s)
% outputs:
%   stability = true/false whether parameterization is stable

r = zeros(length(alpha),1);
for i = 1:length(alpha)
    r(i) = alpha(i)*dt*dx^(-2);
end

if max(r) <= 1/2
    stability = true;
else
    stability = false;
end

end