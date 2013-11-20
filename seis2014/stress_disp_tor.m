function dsdd = stress_disp_tor(r,s)

global omega l rho mu % rspan radius sd

% The input values of s(1) and s(2) are W(r) and T(r) respectively.
% The returned deriatives are stored in dsdd

% structural values at radius r: density and rigidity
% note: use own function here or use globally defined rho and mu (constants)
%[rho,mu] = earthfun(r);

% displacement (first row of equation 1)
dsdd(1,1) = s(1) / r + s(2) / mu;

% stress (second row of equation 2)
dsdd(2,1) = ((l-1)*(l+2)*mu/(r*r) - rho*omega*omega)*s(1) - 3*s(2)/r;
