clear all;
close all;

rng('default');

rho = [1000; 1000; 1000];

m = log(rho);

thk = [20 20];
t = logspace(-6, -3, 31);
nt = length(t);

r = 500.0;

a = 1e-14;

scalefn = @(rho, a) asinh(getVMDLayeredTransient(...
    [t(1) t(end)], r, rho, thk, 0.0, 1 / a));

dHzdt = getVMDLayeredTransient([t(1) t(end)], r, exp(m), thk, 0.0, 1.0);

dm = rand(length(m), 1);
dm = dm ./ norm(dm);

h = logspace(-16, -2, 43);

e0 = zeros(length(h), 1);
e1 = zeros(length(h), 1);

J = getJ('data', dHzdt, 'rho', rho, 'thickness', thk, 'times', t, 'obs', r, ...
    'pert', 1e-8, ...
    'scale_data', @asinh, 'scale_asinh', a, 'protem', false);

for k = 1:length(h)
    dHzdt_p = scalefn(exp(m + h(k) * dm), a);
    e0(k) = norm(dHzdt_p - asinh(dHzdt ./ a));
    e1(k) = norm(dHzdt_p - asinh(dHzdt ./ a) - h(k) * J * dm);
end

%%
figure(1);
loglog(h, e0 ./ e0(end), 'bx', h, h ./ h(end), 'b--', h, e1 ./ e1(end), 'rx', h, h.^2 ./ h(end).^2, 'r--');
grid();
legend('e_0', 'O(h)', 'e_1', 'O(h^2)');
xlabel('h')
ylabel('norm of error functional')
ylim([1e-15 10])
set(gca, 'XDir', 'reverse');
set(gcf, 'position', [1 565 560 420]);
%%
% figure(2);
% imagesc(J);
% colorbar();
% set(gcf, 'position', [1 65 560 420]);