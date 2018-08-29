function [m, obj_fn] = PROTEM_1D_Inversion(varargin)
%PROTEM_1D_Inversion inversion of PROTEM data

p = inputParser;

addParameter(p, 'data', [], @isnumeric);
addParameter(p, 'times', [], @isnumeric);
addParameter(p, 't0', 0.0, @isnumeric);
addParameter(p, 'obs', [], @isnumeric);
addParameter(p, 'depth', [], @isnumeric);
addParameter(p, 'rho', [], @isnumeric);
addParameter(p, 'coil_length',1,@isnumeric);
addParameter(p, 'ng',10,@isnumeric);

addParameter(p, 'thickness', [], @isnumeric);
addParameter(p, 'fix', [], @isnumeric);
addParameter(p, 'verbose', true, @islogical);

addParameter(p, 'maxit', 1e2, @isnumeric);

addParameter(p, 'lambda', 1e2, @isnumeric);
addParameter(p, 'lambda_min', 1e-12, @isnumeric);
addParameter(p, 'reduce_lambda', 1e-1, @isnumeric);
addParameter(p, 'lambda_threshold', 1e-1, @isnumeric);

addParameter(p, 'scale_data', @asinh, @(x) isa(x, 'function_handle'));
addParameter(p, 'scale_asinh', 1e-12, @isnumeric);
addParameter(p, 'goal_normgc', 1e-4, @isnumeric);
addParameter(p, 'goal_obj', 1e-8, @isnumeric);

addParameter(p, 'pert', 1e-6, @isnumeric);

addParameter(p, 'goal_obj_diff', 1e-5, @isnumeric);

addParameter(p, 'linesearch_threshold', 1e-6, @isnumeric);
addParameter(p, 'armijo_c1', 1e-4, @isnumeric);

parse(p, varargin{:});

dobs = p.Results.data;
assert(~isempty(dobs), 'Empty dataset.');

t = p.Results.times;
assert(~isempty(t), 'Empty time range.');

t0 = p.Results.t0;

assert(length(t) == length(dobs), 'Data and times not equal.');

r = p.Results.obs;
assert(r > 0.0, 'r must be > 0.');

z = p.Results.depth;

rho = p.Results.rho;
assert(~isempty(rho), 'Missing starting model layer resistivities.');

thk = p.Results.thickness;
assert(length(thk) == length(rho) - 1, 'Wrong layer thicknesses.');

fix = p.Results.fix;

if isempty(fix)
    fix = zeros(size(rho));
end

assert(length(fix) == length(rho), 'Wrong number of fixed params.');
assert(~all(fix), 'At least one parameter must be free.');

ng = p.Results.ng;
L = p.Results.coil_length;


update = ~fix;
P = eye(length(rho));
if ~all(update)
        P(find(~update), :) = []; %#ok
end
R = eye(min(size(P')));


verbose = p.Results.verbose;
max_gn_iterations = p.Results.maxit;

lambda = p.Results.lambda;
lambda_min = p.Results.lambda_min;
lambda_scaling_factor = p.Results.reduce_lambda;
lambda_droptol = p.Results.lambda_threshold;

scale_data = p.Results.scale_data;
a = p.Results.scale_asinh;
pert = p.Results.pert;

scalefn = @(x, a) scale_data(x ./ a);

c1 = p.Results.armijo_c1;
ls_drop = p.Results.linesearch_threshold;

tol = p.Results.goal_normgc;
obj_fn_goal = p.Results.goal_obj;
goal_obj_diff = p.Results.goal_obj_diff;

rho_ref = 1 * ones(size(rho));

% Find time interval boundaries
%
[tmin, tmax] = minmax(t);
lmi = floor(log10(tmin));
lma = ceil(log10(tmax));
nt_ = 1 + 10 * (lma - lmi);
tstart = 10^lmi;
tend = 10^lma;

tref = logspace(floor(log10(tmin)), ceil(log10(tmax)), nt_);


% Some useful function handles
%



[x,w] = gauleg(-L/2,L/2,ng);
getData = @(rho) simulate_PROTEMSQLOOP(t, r, rho, thk, t0, z,x,w,L);

residual = @(dobs, d, a) scalefn(dobs, a) - scalefn(d, a);
data_norm = @(r) 0.5 * norm(r)^2;
model_norm = @(W, m, lambda) 0.5 * lambda * m' * W' * W * m;
objective_function = @(r, W, m, lambda) data_norm(r) + 0 * model_norm(W, m, lambda);

% First (G, gradient) and second (W) derivative of model parameters wrt to
% depth, as well as identity matrix (I) of compatible dimension
%
W = eye(length(rho));

% Transform model parameters, here we take the natural log of rho
%
log_p = log(rho);
log_p_ref = log(rho_ref);

% Calculate response of starting model
%
d = getData(rho);

% Calculate data residual for starting model
%
data_residual = residual(dobs, d, a);

obj_fn_current = ...
    data_norm(data_residual) + ...
    model_norm(W, log_p - log_p_ref, lambda);

% Initialize GN monitoring
%
obj_fn = zeros(max_gn_iterations, 1);
norm_gc = Inf();
diff_obj_fn = Inf();

gn_it_counter = 0;

obj_fn(1) = obj_fn_current;

% GN starts here
%

while ( ...
        norm_gc) > tol && ...
        (gn_it_counter < max_gn_iterations) && ...
        (obj_fn_current > obj_fn_goal) && ...
        (diff_obj_fn > goal_obj_diff)
    %
    % We have data d, current model log_p, and have evaluated the objective
    % function.
    %
    % Jacobian J = du/dm * dm / drho
    J = getJ('data', d, ...
        'rho', exp(log_p), ...
        'thickness', thk, ...
        'times', t, ...
        'obs', r, ...
        'depth',z,...
        'scale_asinh', a, ...
        'pert', pert, ...
        'protem', true, ...
        't0', t0);
    
    J = J * P';
    
    % Gradient
    gc = J' * data_residual;
    norm_gc = norm(gc);
    
    % GN search direction
    %
    gn_dir = (J' * J + lambda * R) \ gc;
    gn_dir = P' * gn_dir;
    
    if verbose
        fprintf('# %3d | ', gn_it_counter);
        fprintf('JTr=%1.3e | ', norm_gc);
    end
    
    % LS
    step_size = 1.0;
    % take a fraction of the GN direction and evaluate the response d_
    log_p_ = log_p + step_size * gn_dir;
    d_ = getData(exp(log_p_));
    
    obj_fn_tmp = objective_function(residual(dobs, d_, a), W, log_p_ - log_p_ref, lambda);
    diff_obj_fn = obj_fn_current - obj_fn_tmp;
    
    % Directional derivative for Armijo rule
    dir_deriv = gn_dir' * (P' * gc);
    
    gn_it_counter = gn_it_counter + 1;
    
    while (diff_obj_fn < -c1 * step_size * dir_deriv)
        step_size = step_size / 2.0;
        if step_size < ls_drop
%             error('Line search breakdown. Exit.');
                break;
        end
        log_p_ = log_p + step_size * gn_dir;
        d_ = getData(exp(log_p_));
        obj_fn_tmp = objective_function(residual(dobs, d_, a), W, log_p_ - log_p_ref, lambda);
        diff_obj_fn = obj_fn_current - obj_fn_tmp;
    end
    
    % Update
    rho = exp(log_p_);
    log_p = log_p_;
    d = d_;
    data_residual = residual(dobs, d, a); %scalefn(dobs, a) - scalefn(d, a);
    obj_fn_current = obj_fn_tmp;
    obj_fn(gn_it_counter) = obj_fn_current;
    
    % Reduce lambda when difference in objective function falls below droptol
    if abs(diff_obj_fn) < lambda_droptol
        lambda = lambda * lambda_scaling_factor;
    end
    
    if lambda < lambda_min
        lambda = lambda_min;
    end
    if verbose
        fprintf('alpha=%1.2e | lambda=%1.2e | obj_fn=%1.3e | diff_obj=%1.3e\n', ...
            step_size, lambda, obj_fn_current, diff_obj_fn);
        
        figure(1);
        loglog(t, abs(dobs), t, abs(d));
        plotTransient(t, dobs, [], [1 0 0]);
        hold on
        plotTransient(t, d, [], [0 0 1]);
        %     grid('on');
        legend('dobs', '', 'd', '');
        hold off
        set(gcf, 'position', [1 565 560 420]);
        figure(2);
        plot_model(rho, thk);
        set(gcf, 'Position', [613 549 560 419]);
        drawnow();
    end
end

% Return value
m = rho;
obj_fn = obj_fn(1:gn_it_counter);
