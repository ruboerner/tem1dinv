close all;
clear all;
clc();

rng('default');

t = logspace(-6, -2, 41);
t = t(:);

rho = [100; 10; 100; 1000];
thk = [50; 10; 30];
r = 100.0;

dobs = getVMDLayeredTransient([t(1) t(end)], r, rho, thk, 0.0, 1.0);
dobs = dobs + t.^(0.5)  .* randn(length(dobs), 1) * 1e-14;

nl = 10;
rho_guess = 100 * ones(nl, 1);
% rho_guess(1) = 100;
thk_guess = 10 * ones(nl - 1, 1);
fixp = zeros(nl, 1);
fixp([1]) = 1;

%%
ti = tic();

[rho_new, obj_fn] = TEM_1D_Inversion( ...
    'data', dobs, ...
    'obs', r, ...
    'times', t, ...
    'rho', rho_guess, ...
    'thickness', thk_guess, ...
    'fix', fixp, ...
    'scale_asinh', 1e-12, ...
    'goal_obj', 1e-5, ...
    'goal_obj_diff', 3e-7, ...
    'lambda_min', 1e-10, ...
    'pert', 1e-6, ...
    'verbose', true);

tInv = toc(ti);

fprintf('Inversion completed after %.2f s and %d iterations.\n', tInv, length(obj_fn));
%%
figure(2);
plot_model(rho, thk);
hold on
plot_model(rho_new, thk_guess);
hold off
legend('true', 'est');
xlim([1 1000]);
set(gcf, 'Position', [613 549 560 420]);

%%
figure(3);
semilogy(obj_fn ./ obj_fn(1));
xlabel('Iterations');
ylabel('Rel. objective function');
grid();
xlim([1 length(obj_fn) + 1]);
set(gcf, 'position', [52 43 550 420]);
% matlab2tikz('obj_fn.tex');
