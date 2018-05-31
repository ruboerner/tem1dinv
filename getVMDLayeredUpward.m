function b = getVMDLayeredUpward(u, f, typ, rho, d)
% TODO: Update help text.
%
% FUNCTION B=BTP(U,F,TYP,RHO,D)
% Admittance of a layered halfspace
% for TE (TYP=1) and TM (TYP=2) mode.
% U: wave number (1/m)
% F: frequency (1/s)
% RHO: layer resitivities
% D: layer thicknesses

% make sure resistivities and layer thicknesses are 3-vectors
nl = length(rho);
rho = shiftdim(rho(:), -2);
d = shiftdim(d(:), -2);

% constants
mu0 = 4e-7 * pi;
c = 1i * mu0 * 2 * pi * f;

% homogeneous half-space
b = sqrt(bsxfun(@plus, u .^ 2, c / rho(nl)));

% layered half-space
if nl > 1
    beta = 1;
    for nn = nl - 1:-1:1
        alpha = sqrt(bsxfun(@plus, u .^ 2, c / rho(nn)));
        if typ == 2
            beta = rho(nn) / rho(nn + 1);
        end
        cth = exp(-2 * d(nn) * alpha);
        cth = (1 - cth) ./ (1 + cth);
        b = (b + alpha .* beta .* cth) ./ (beta + cth .* b ./ alpha);
    end
end
