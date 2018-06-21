function dBzdt = simulate_PROTEMSQLOOP(t, r, rho, thk, t0, z,x,w,L)
%simulatePROTEM(t, r, rho, thk, t0)
[tmin, tmax] = minmax(t);
lmi = floor(log10(tmin));
lma = ceil(log10(tmax));
nt_ = 1 + 10 * (lma - lmi);
tstart = 10^lmi;
tend = 10^lma;
tref = logspace(floor(log10(tmin)), ceil(log10(tmax)), nt_);
dBzdt = zeros(nt_,1);


for i = 1:length(x)
    for j = 1:length(x)
        xx = x(i)+r;
        yy = x(j);
        dBzdtgauss = getVMDLayeredTransient([tstart,tend],r,rho,thk,z);
        dBzdt = dBzdt+w(i)*w(j)*dBzdtgauss;
    end
end
dBzdt = dBzdt/(L*L);

v = interp_transient( ...
        tref, ...
        dBzdt, ...
        t);
dBzdt = correctRampTime(t, v, t0);
