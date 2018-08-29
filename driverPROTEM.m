clear
close all
clc
warning off
%%
%for Tharandt
load('L:\TU-Freiberg\TEM\Reproduzierbarkeit\ThaWa\ALL DATA\transients3_stacked.mat')
L = 100;
r = 270;
%{
% for Canada
load('L:\TU-Freiberg\TEM\2014_Alberta\TEM\Messungen TEM 07_2016\DATA_Reprocessed\Alberta_stacked.mat')
load('L:\TU-Freiberg\TEM\2014_Alberta\TEM\Messungen TEM 07_2016\DATA_Reprocessed\GPSutm.mat');
records = data_stacked(ismember(data_stacked.Loop, 'NE'), :);
records = records(ismember(records.transmitter,57),:);
coil = cell2mat([GPS.Lon(31:34), GPS.Lat(31:34)]);
obs  = cell2mat([GPS.Lon(2), GPS.Lat(2)]);
% check distace from midpoint of loop to obs
MM = mean(coil);
r = norm(-MM+obs);
L = 400;
%}
%% data handling
%for Tharandt
id = 22;
range = 1:25;
records = data_stacked(id,:);
dobs = records.Z_field{1}(range);
z = 0;
t = records.t{1}(range);



% for Canada 
%{
% id in record table
id = 22;
% specify gates
range = 7:30;
% truncate transient according to range
t    = records.t{id}(range);
dobs = -records.Z_field{id}(range);
dobs_err = records.Z_err{id}(range);
% depth of observation
z    = records.depth(id); 
% z = 275;
%}

%% model specification
% thk_guess = 2 * ones(length(0:5:z), 1);
% thk_guess = [thk_guess; logspace(1,log10(250),15)']
% rho_guess = 10 * ones(length(thk_guess)+1, 1);
% thk_guess = 50 * ones(nl - 1, 1);
nl = 30;
rho_guess = 100 * ones(nl, 1);
% thk_guess = logspace(1,log10(250),nl-1);
thk_guess = 40 * ones(nl-1,1);
% cumsum(thk_guess)
% fix = zeros(length(thk_guess)+1,1);
fix = zeros(nl,1);
% fix(1) = 1;
% rho_guess(1) = 1e3;
%%
ti = tic();
[rho_new, obj_fn] = PROTEM_1D_Inversion( ...
    'coil_length',L,...
    'data', dobs, ...
    'ng',3,...
    't0', 425e-6, ...
    'obs', r, ...
    'depth',z, ...
    'times', t, ...
    'rho', rho_guess, ...
    'thickness', thk_guess, ...
    'fix', fix, ...
    'scale_asinh', 1e-16, ...
    'goal_obj', 1e-4, ...
    'goal_obj_diff', 1e-4, ...
    'lambda_threshold', 1e-2, ...
    'lambda_min',1e-8,...
    'verbose', true,...
    'reduce_lambda',0.1,...
    'maxit',100);
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
