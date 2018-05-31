function vc = correctRampTime(t, v, t0)
%function vc = correctRampTime(t, v, t0)
%
%  Corrects a given impulse response v(t) for a finite shut-off
%  ramp time t0 (for details, cf. Fitterman, Anderson (1987)).
%
%  Basically, the function computes the integral
%
%       vc(t) = 1 / t0 \int\limits_{t}^{t + t0} v(t') dt'
%
%  for any t > 0 with the help of an adaptive Gauss quadrature technique.
%
%  It is assumed that the linear ramp ends at t=0, i.e., the ramp starts at
%  t = -t0.
%
%  RUB (2013)
%
pp = spline(t, v);

interpolant = @(tt) ppval(pp, tt);

vc = zeros(size(v));
for kk = 1:length(t)
    vc(kk) = integral(interpolant, ...
        t(kk),t(kk) + t0) / t0;
end
