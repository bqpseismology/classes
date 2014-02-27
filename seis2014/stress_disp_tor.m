function dWT = stress_disp_tor(r,WT)

global omega l rho mu % rspan radius WT

% The input values of WT(1) and WT(2) are W(r) and T(r) respectively.
% The returned deriatives are stored in dWT

% structural values at radius r: density and rigidity
% note: use own function here or use globally defined rho and mu (constants)
%[rho,mu] = earthfun(r);

% displacement (first row of equation 1)
dWT(1,1) = WT(1) / r + WT(2) / mu;

% stress (second row of equation 2)
dWT(2,1) = ((l-1)*(l+2)*mu/(r*r) - rho*omega*omega)*WT(1) - 3*WT(2)/r;
