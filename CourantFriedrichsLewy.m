function [stability] = CourantFriedrichsLewy(E, dt, dx)
%% evaluate Courant-Friedrichs-Lewy 1D stability criterion
% inputs:
%   E = vector of bank erosion rates (m/s)
%   dt = timestep (s)
%   dx = spatial step (m)
% output:
%   stability = t/f if model parameterization is stable

C = zeros(length(E),1);
for i = 1:length(C)
    C(i) = E(i)*dt/dx;
end

if max(C) < 1
    stability = true;
else
    stability = false;
end

end