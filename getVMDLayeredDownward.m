function [a, ap, b1] = getVMDLayeredDownward(u, f, rho, d, z)
% TODO: Update help text.
%
% [a,ap,b1]=downward(u,f,rho,d,z)
%
% Liefert Korrekturfaktoren, mit denen das elektromagnetische
% Feld im Punkt z im Innern des Leiters (geschichteter Halbraum)
% auf das Oberflaechenfeld in z=0 bezogen werden kann.
%
% Input:
% u:   Raeumliche Wellenzahl in 1/Meter
% f:   Frequenz in Hertz
% rho: Spez. el. Widerstaende der Schichten, rho(nl) in
%      Ohm*Meter
% d:   Maechtigkeiten der Schichten, d(nl-1) in Meter)
% z:   Punkt im Innern des Leiters in Meter
%
% Output:
% a:   Verhaeltnis der Partialwellenamplitude A(z,u)/A(0,u) (z.B. fuer
%      E_phi, H_z)
% ap:  Verhaeltnis der Partialwellenamplitude A'(z,u)/A'(0,u) (z.B.
%      fuer H_r)
% b1:  Admittanz an der Oberflaeche des geschichteten Halbraums
%
% Ralph-Uwe Boerner (2005)

% make sure resistivities and layer thicknesses are a 3-vectors
nl = length(rho);
rho = shiftdim(rho(:), -2);
d = shiftdim(d(:), -2);

% constants
mu0 = 4e-7 * pi;
iwm = 1i * 2 * pi * f * mu0;

% propagation constant
alpha = sqrt(bsxfun(@plus, u .^ 2, bsxfun(@rdivide, iwm, rho)));

% Oberflächenadmittanz
if nl == 1
    % homogeneous half-space
    b1 = alpha;
    [a, ap] = deal(exp(-alpha .* z));
elseif nl > 1
    % layered half-space
    talphad = tanh(bsxfun(@times, alpha(:, :, 1:nl - 1), d));
    
    % Rekursive Bestimmung der Admittanzen an der Oberkankte der Schicht
    b = complex(zeros(size(alpha)));
    b(:, :, nl) = alpha(:, :, nl);
    for nn = nl - 1:-1:1
        b(:, :, nn) = ...
            alpha(:, :, nn) .* (b(:, :, nn + 1) + alpha(:, :, nn) .* talphad(:, :, nn)) ./ ...
            (alpha(:, :, nn) + b(:, :, nn + 1) .* talphad(:, :, nn));
    end
    
    % reciprocal values gives impedance
    c = 1 ./ b;
    
    % surface admittance for layered half-space
    b1 = b(:, :, 1);
    
    % variation from layer boundary to layer boundary
    [aa, aap] = deal(complex(zeros(size(talphad))));
    for nn = 1:nl - 1
        aa(:, :, nn) = ...
            (b(:, :, nn) + alpha(:, :, nn)) ./ (b(:, :, nn + 1) + alpha(:, :, nn)) .* ...
            exp(-alpha(:, :, nn) * d(nn));
        aap(:, :, nn) = ...
            (1 + alpha(:, :, nn) .* c(:, :, nn)) ./ (1 + alpha(:, :, nn) .* c(:, :, nn + 1)) .* ...
            exp(-alpha(:, :, nn) * d(nn));
    end
    
    % location of layer boundaries
    h = cat(3, 0, cumsum(d));
    
    % find layer index for given z
    ind = find(z >= h(1:end - 1) & z < h(2:end), 1, 'first');

    if ~isempty(ind) % means: ind < nl
        % a layer of finite thickness
        a = 0.5 * prod(aa(:, :, 1:ind - 1), 3) .* (1 + b(:, :, ind) ./ alpha(:, :, ind)) .* (...
            exp(-alpha(:, :, ind) * (z - h(ind))) - ...
            (b(:, :, ind + 1) - alpha(:, :, ind)) ./ (b(:, :, ind + 1) + alpha(:, :, ind)) .* ...
            exp(-alpha(:, :, ind) * (d(ind) + h(ind + 1) - z)));
        ap = 0.5 * prod(aap(:, :, 1:ind - 1), 3) .* (1 + alpha(:, :, ind) .* c(:, :, ind)) .* (...
            exp(-alpha(:, :, ind) * (z - h(ind))) + ...
            (1 - alpha(:, :, ind) .* c(:, :, ind + 1)) ./ (1 + alpha(:, :, ind) .* c(:, :, ind + 1)) .* ...
            exp(-alpha(:, :, ind) * (d(ind) + h(ind + 1) - z)));
    else
        % lower half-space
        ind = nl;
        a = prod(aa, 3) .* exp(-alpha(:, :, ind) .* (z - h(ind)));
        ap = prod(aap, 3) .* exp(-alpha(:, :, ind) .* (z - h(ind)));
    end
else
    error('NL < 1!');
end
