function [Cf, ustar, Hsk] = ManningStrickler(U, S, ks)

g = 9.81;                       % gravitational acceleration (m2/s)
Hsk = (U/8.32)^(2/3) * ks^(1/4) * (g*S)^(-3/4);
ustar = sqrt(g*Hsk*S);
Cf = (ustar / U)^2;

end