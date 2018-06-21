function dBzdt = simulate_PROTEM(t, r, rho, thk, t0, z)
%simulatePROTEM(t, r, rho, thk, t0)
[tmin, tmax] = minmax(t);
lmi = floor(log10(tmin));
lma = ceil(log10(tmax));
nt_ = 1 + 10 * (lma - lmi);
tstart = 10^lmi;
tend = 10^lma;
tref = logspace(floor(log10(tmin)), ceil(log10(tmax)), nt_);
v = interp_transient( ...
        tref, ...
        getVMDLayeredTransient([tstart tend], r, rho, thk, z, 1.0), ...
        t);
dBzdt = correctRampTime(t, v, t0);
