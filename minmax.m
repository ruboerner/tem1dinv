function [mi, ma] = minmax(v)
%minmax Smallest and largest element.
%
% [MI, MA] = minmax(V) computes the smallest and largest values of V.
%
mi = min(v(:));
ma = max(v(:));
if nargout == 0
    fprintf('%.6e ... %.6e\n', mi, ma);
end