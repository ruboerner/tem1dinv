clear all
close all
clc()

%%
load data/transients_stacked.mat

records = data_stacked(ismember(data_stacked.date, '21.02.'), :);

t    = records.t{1};
dobs = records.Z_field{12};
r    = 280; % FIXME

% t = t(10:end);
% dobs = dobs(10:end);

%%
nl = 12;
rho_guess = 100 * ones(nl, 1);
thk_guess = 10 * ones(nl - 1, 1);

%%
ti = tic();
[rho_new, obj_fn] = PROTEM_1D_Inversion( ...
    'data', dobs, ...
    't0', 57e-6, ...
    'obs', r, ...
    'times', t, ...
    'rho', rho_guess, ...
    'thickness', thk_guess, ...
    'fix', [1 0 0 0 0 0 0 0 0 0 0 0], ...
    'scale_asinh', 1e-12, ...
    'goal_obj', 1e-2, ...
    'goal_obj_diff', 1e-4, ...
    'lambda_threshold', 1e-2, ...
    'verbose', true);
tInv = toc(ti);
fprintf('Inversion completed after %.2f s and %d iterations.\n', tInv, length(obj_fn));

%%
figure(3);
semilogy(obj_fn ./ obj_fn(1));
xlabel('Iteration');
ylabel('Rel. objective function');
grid();
xlim([1 length(obj_fn) + 1]);
set(gcf, 'position', [52 43 550 420]);
