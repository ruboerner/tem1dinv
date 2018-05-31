function d = interp1_transient(t, v, t_)
%interp1_transient(t, v, t_)
%
F = griddedInterpolant(t, v);
d = F(t_);