function varargout = getVMDLayeredHarmonic(outs, f, r, rho, d, z, h, dipm)
% TODO: Update help text.
%
% [FIELDS,RR]=VMD_F(IOUT,F,RHO,D,Z,H,RMIN,NR,TM)
%
% Berechnung der Felder eines zeitharmonischen Vertikalen
% Magnetischen Dipols (VMD) ueber einem geschichteten Halbraum
% mittels Hankeltransformation.
% Verwendet werden Zylinderkoordinaten mit z positiv nach
% unten.
% Der VMD befindet sich in r=0 und z=-h, h>0. Die
% Horizontalkomponenten der Aufpunkte RR sind logarithmisch
% aequidistant verteilt und werden aus RMIN und NR berechnet. Die
% Felder werden fuer alle z bestimmt.
% Die verwendeten Filterkoeffizienten sind fuer 10 Stuetzstellen
% pro Dekade vorgegeben.
% Der geschichtete Halbraum wird durch Angabe von NL
% spezifischen Widerstaenden und NL-1 Maechtigkeiten beschrieben.
%
% Input:
% IOUT:  Auswahl der zu berechnenden Felder
%        IOUT(3)
%        Vorschrift: [HZ,HR,EPHI].*IOUT
%        Beispiele: IOUT=[1,0,0]: Berechnung von HZ
%                   IOUT=[0,1,1]: Berechnung von HR und EPHI
% F:     Frequenz in Hertz
% RHO:   Spezifische Widerstaende in Ohm*m, RHO(NL)
% D:     Maechtigkeiten in m, D(NL-1), [] fuer hom. Halbraum
% Z:    z-Koordinate des Empfaengers in Meter (default: ZE=0)
% H:    Hoehe des Senders in Meter (default: H=0), H>=0
% RMIN:  Kleinster Empfaengerabstand in Meter (default: RMIN=1)
% NR:    Maximale Anzahl von Entfernungen (default: NR=41)
% TM:    Magnetische Dipolmoment in Ampere*Meter^2 (default:
%        TM=1)
%
% Output:
% In Abhaengigkeit von IOUT:
% FIELDS(1,NR): HZ
% FIELDS(2,NR): HR
% FIELDS(3,NR): EPHI
% HZ:    Vertikalkomponente des Magnetfeldes in Ampere/Meter,
%        Typ COMPLEX
% HR:    Radialkomponente des Magnetfeldes in Ampere/Meter,
%        Typ COMPLEX
% EPHI:  Tangentialkomponente des elektrischen Feldes in
%        Volt/Meter, Typ COMPLEX
% RR:    Aus RMIN und NRMAX gebildete Entfernungen RR(NR)
%
% Beispiel:
%  [FIELDS,RR]=VMD_F([1,0,0],1000,300,[],0,0,1,21,1);
%  berechnet HZ-Komponente fuer F=1000 Hertz, 300 Ohm*m Halbraum,
%  Sender in z=0, Empfaenger in z=0, RMIN=1 Meter, 21 Abstaende,
%  Magnetisches Dipolmoment 1 Ampere*Meter^2
%
% Aufgerufene Funktionen: HANKELFC, BTP, DOWNWARD
%
% R.-U. Boerner
%

% check parameters
assert(nargin >= 6, ...
    'Expected at least six parameters.');

assert(iscellstr(outs), ...
    'Expected cell array of strings naming the requested outputs.');

assert(length(outs) == nargout, ...
    'Expected number of outputs to match the number of requested outputs.');

% provide defaults
if nargin < 7
    h = 0;
end

if nargin < 8
    dipm = 1;
end

% check and handle radii
if isscalar(r)
    rmin = r;
    nr = 1;
elseif length(r) == 2
    assert(r(1) < r(2), ...
        'Expected range for ''r'' to be monotonic.');
    rmin = 10 ^ floor(log10(r(1)));
    rmax = 10 ^ ceil(log10(r(2)));
    nr = 10 * log10(rmax / rmin) + 1;
else
    error('Expected parameter ''r'' to be a scalar or two-element vector.');
end

% check position of transmitter
assert(h >= 0, ...
    'Transmitter only allowed to be at z <= 0.');

% constants
mu0 = 4e-7 * pi;

% radii (row vectors)
q = 10 ^ 0.1;
rr = rmin * q .^ ((1:nr) - 1);
R = sqrt(rr .^ 2 + (z + h) .^ 2);

% filter coefficients for wave numbers/radii (row vectors)
[fc0, nc, nc0] = getHankelFC('j0');
fc1 = getHankelFC('j1');
[fc0, fc1] = deal(fc0.', fc1.');
ncnr = nc + nr - 1;

% wave numbers (row vector)
nu = 1:ncnr;
n = nc0 - nc + nu;
q = 0.1 * log(10);
u = exp(-(n - 1) * q) / rmin;

% frequencies (column vector)
nf = length(f);
f = f(:);

% admittances for all wave numbers and frequencies
if z <= 0
    % air half-space
    m1 = getVMDLayeredUpward(u, f, 1, rho, d);
else
    % conductive half-space
    [aa, aap, m1] = getVMDLayeredDownward(u, f, rho, d, z);
end

% Kernfunktionen
if z <= 0
    % air half-space
    e = exp(u .* (z - h));
    e = repmat(e, [nf, 1]);
    u = repmat(u, [nf, 1]);
    delta = (u - m1) ./ (u + m1) .* e;
else
    % conductive half-space
    e = exp(-u .* h);
    e = repmat(e, [nf, 1]);
    u = repmat(u, [nf, 1]);
    delta = 2 * u .^ 2 ./ (u + m1) .* e;
end

% signature of fields in wave number domain
if z <= 0
    % air half-space
    sig_Ef = delta .* u;
    sig_Hz = sig_Ef .* u; % = delta .* u .* u;
%     sig_Hr = -sig_Hz; % = -delta .* u .* u;
else
    % conductive half-space
    sig_Ef = delta .* aa;
    sig_Hz = sig_Ef .* u; % = delta .* u .* aa;
%     sig_Hr = delta .* m1 .* aap;
end

% convolution (wave numbers -> radii)
rr_inv = repmat(1 ./ rr, [nf, 1]);
Hz = rr_inv .* conv2(1, fc0, sig_Hz, 'valid');
% Hr = rr_inv .* conv2(1, fc1, sig_Hr, 'valid');
% Ef = rr_inv .* conv2(1, fc1, sig_Ef, 'valid');

% absolute field values for air-half space
if z <= 0
    % field values in air
    air_Hz = (3 * (z + h) .^ 2 - R .^ 2) ./ R .^ 5;
%     air_Hr = 3 * rr * (z + h) ./ R .^ 5;
%     air_Ef = rr ./ R .^ 3;

    % add field values in air (not normalized)
    Hz = bsxfun(@plus, Hz, air_Hz);
%     Hr = bsxfun(@plus, Hr, air_Hr);
%     Ef = bsxfun(@plus, Ef, air_Ef);
end

% auxiliary quantities for normalization
bh = 1 / (4 * pi); % A/m oder nT
% omega = 2 * pi * f;
% cc = 1i * omega * mu0;

% normalization of field values
Hz = (dipm * bh) * Hz;
% Hr = (dipm * bh) * Hr;
% Ef = (-dipm / (4 * pi)) * bsxfun(@times, cc, Ef);

% construct output
varargout = cell(1, length(outs));
for i = 1:length(outs)
    switch outs{i}
        case 'r'
            varargout{i} = rr;
        case 'Hz'
            varargout{i} = Hz;
        case 'Hr'
            varargout{i} = Hr;
        case 'Ef'
            varargout{i} = Ef;
        otherwise
            error('Unrecognized requested output ''%s''.', outs{i});
    end
end
