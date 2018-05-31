function dHzdt = getVMDLayeredTransient(t, r, rho, d, z, dipm)


% provide defaults
if nargin < 6
    dipm = 1;
end

% check and handle times
if isscalar(t)
    tmin = t;
    nt = 1;
elseif length(t) == 2
    assert(t(1) < t(2), ...
        'Expected range for ''t'' to be monotonic.');
    tmin = 10 ^ floor(log10(t(1)));
    tmax = 10 ^ ceil(log10(t(2)));
    nt = 10 * log10(tmax / tmin) + 1;
else
    error('Expected parameter ''t'' to be a scalar or two-element vector.');
end

% parameters r, rho, d, and z will be checked in getVMDLayeredHarmonic

% times (column vector)
q = 10 ^ 0.1;
tt = tmin * q .^ ((1:nt).' - 1);

% filter coefficients for frequencies/times (column vector)
[fcs, nc, nc0] = getHankelFC('cos');
ncnt = nc + nt - 1;

% angular frequencies (column vector)
omega = q .^ (1 - (-nc + nc0 + (1:ncnt).')) ./ tmin;
omega_sqrt = sqrt(omega);

% array sizes:
% * frequency domain quantities: [ncnt, nr]
% * time domain quantities: [nt, nr]

% compute raw fields in frequency domain
raw_Hz = getVMDLayeredHarmonic(...
    {'Hz'}, ...
    omega ./ (2 * pi), r, rho, d, z);

% frequency domain: time derivatives of fields
dHz_f = real(bsxfun(@times, raw_Hz, omega_sqrt));
% dHr_f = real(bsxfun(@times, raw_Hr, omega_sqrt));
% dEf_f = real(bsxfun(@times, raw_Ef, omega_sqrt));

% frequency domain: fields
% Hz_f = -imag(bsxfun(@times, raw_Hz, 1 ./ omega_sqrt));
% Hr_f = -imag(bsxfun(@times, raw_Hr, 1 ./ omega_sqrt));
% Ef_f = -imag(bsxfun(@times, raw_Ef, 1 ./ omega_sqrt));

% convolution: time derivatives of fields (frequencies -> times)
dHz = conv2(fcs, 1, dHz_f, 'valid');
% dHr = conv2(fcs, 1, dHr_f, 'valid');
% dEf = conv2(fcs, 1, dEf_f, 'valid');

% convolution: fields (frequencies -> times)
% Hz = conv2(fcs, 1, Hz_f, 'valid');
% Hr = conv2(fcs, 1, Hr_f, 'valid');
% Ef = conv2(fcs, 1, Ef_f, 'valid');

% scale results: compute adjustment factor
fac = dipm * sqrt(2 ./ (pi * tt));

% scale results: adjust time derivatives of fields
dHz = bsxfun(@times, fac, dHz);
% dHr = bsxfun(@times, fac, dHr);
% dEf = bsxfun(@times, fac, dEf);

% scale results: adjust fields
% Hz = bsxfun(@times, fac, Hz);
% Hr = bsxfun(@times, fac, Hr);
% Ef = bsxfun(@times, fac, Ef);

dHzdt = dHz * pi * 4e-7;
