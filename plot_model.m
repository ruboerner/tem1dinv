function plot_model(rho, thk, marker)

if nargin < 3
    marker = '.-';
end
r = reshape(repmat(reshape(rho, 1, []), 2, 1), [], 1);
z = [0; reshape(repmat(reshape(cumsum(thk), 1, []), 2, 1), [], 1); 1e3];

plot(r, z, marker);
ylim([0 1.4 * z(end - 1)]);
xlabel('rho in Ohmm');
ylabel('depth in m');
grid('on');
set(gca, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'left', ...
    'YDir', 'reverse', ...
    'XScale', 'log');