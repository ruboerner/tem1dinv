function J = getJ(varargin)
%getJ Jacobian dBzdt'(rho) for the 1D TEM case
%

p = inputParser;

addParameter(p, 'data', [], @isnumeric);
addParameter(p, 'times', [], @isnumeric);
addParameter(p, 't0', 0.0, @isnumeric);
addParameter(p, 'obs', [], @isnumeric);
addParameter(p, 'depth', 0, @isnumeric);
addParameter(p, 'rho', [], @isnumeric);
addParameter(p, 'thickness', [], @isnumeric);
addParameter(p, 'scale_data', @asinh, @(x) isa(x, 'function_handle'));
addParameter(p, 'scale_asinh', 1e-12, @isnumeric);
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), 'Parameter must be > 0');
addParameter(p, 'pert', 1e-6, validationFcn);
addParameter(p, 'protem', true, @islogical);

parse(p, varargin{:});

protem = p.Results.protem;
a = p.Results.scale_asinh;
rho = p.Results.rho;
thk = p.Results.thickness;
r = p.Results.obs;
z = p.Results.depth;
t = p.Results.times;
d = p.Results.data;
t0 = p.Results.t0;
scale_data = p.Results.scale_data;
pert = p.Results.pert;
scalefn = @(x, a) scale_data(x ./ a);

if protem
    getData = @(rho) simulate_PROTEM(t, r, rho, thk, t0, z);
else
    tstart = t(1);
    tend = t(end);
    getData = @(rho) getVMDLayeredTransient([tstart tend], r, rho, thk, 0.0, 1.0);
end

d = scalefn(d, a);
J = zeros(length(t), length(rho));

for i = 1:length(rho)
    pr = rho;
    pr(i) = (1 + pert) * rho(i);
    pdH = getData(pr);
    pdH = scalefn(pdH, a);
    J(:, i) = (pdH - d) / (pert * rho(i)) * rho(i);
end