function plot_model(rho, thk)

r = reshape(repmat(reshape(rho, 1, []), 2, 1), [], 1);
lb = cumsum(thk); %find lower boundary of model and extend
z = [0; reshape(repmat(reshape(cumsum(thk), 1, []), 2, 1), [], 1); 2*lb(end)];

plot(r, z, '.-');
ylim([0 1.4 * z(end - 1)]);
xlabel('rho in Ohmm');
ylabel('depth in m');
grid('on');
set(gca, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'left', ...
    'YDir', 'reverse', ...
    'XScale', 'log');