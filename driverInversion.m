close all;
clearvars;
clc();

rng('default');

t = logspace(-6, -2, 41);
t = t(:);

rho = [100; 10; 100];
thk = [50; 30];
r = 100.0;

dobs = getVMDLayeredTransient([t(1) t(end)], r, rho, thk, 0.0, 1.0);
% dobs = dobs + t.^(-0.5)  .* randn(length(dobs), 1) * 1e-15;

nl = 12;
rho_guess = 500 * ones(nl, 1);
% rho_guess(8) = 100;
rho_guess(1) = 100;
thk_guess = 10 * ones(nl - 1, 1);
fixp = zeros(nl, 1);
% fixp([8]) = 1;

%%
ti = tic();

[rho_new, obj_fn, lambda_history] = TEM_1D_Inversion( ...
    'data', dobs, ...
    'obs', r, ...
    'times', t, ...
    'rho', rho_guess, ...
    'thickness', thk_guess, ...
    'fix', fixp, ...
    'scale_asinh', 1e-12, ...
    'goal_obj', 1e-6, ...
    'goal_obj_diff', 1e-8, ...
    'lambda', 1e0, ...
    'reduce_lambda', 1e-1, ...
    'lambda_min', 1e-8, ...
    'lambda_threshold', 0.5, ...
    'armijo_c1', 1e-4, ...
    'pert', 1e-6, ...
    'verbose', true);

tInv = toc(ti);

fprintf('Inversion completed after %.2f s and %d iterations.\n', tInv, length(obj_fn));
%%
figure(2);
plot_model(rho, thk, '.-');
hold on
plot_model(rho_new, thk_guess, '-');
hold off
legend('true', 'est');
xlim([1 1000]);
set(gcf, 'Position', [613 50 560 420]);

%%
figure(3);
it_number = 0:(length(obj_fn) - 1);
yyaxis left
% semilogy(obj_fn ./ obj_fn(1));
semilogy(it_number, obj_fn, '.-', 'MarkerSize', 12);
axl = gca();
xlabel('Iteration number');
ylabel('Rel. objective function');
% grid();
set(axl,'YGrid', 'on', 'YMinorGrid', 'on')
set(axl,'XGrid', 'on', 'XMinorGrid', 'on')

xlim([0 length(obj_fn)]);
yyaxis right
semilogy(it_number, lambda_history, '.-', 'MarkerSize', 12);
% grid();
ylabel('Regularization paramater');
set(gcf, 'position', [1191 50 550 420]);
% matlab2tikz('obj_fn.tex');
