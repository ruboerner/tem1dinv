function [m, obj_fn, lambda_history] = TEM_1D_Inversion(varargin)
%TEM_1D_Inversion inversion of synthetic TEM data
%

p = inputParser;

addParameter(p, 'data', [], @isnumeric);
addParameter(p, 'times', [], @isnumeric);
addParameter(p, 'obs', [], @isnumeric);
addParameter(p, 'rho', [], @isnumeric);

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
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0), 'Parameter must be > 0');
addParameter(p, 'pert', 1e-6, validationFcn);
addParameter(p, 'goal_obj_diff', 1e-5, @isnumeric);

addParameter(p, 'linesearch_threshold', 1e-6, @isnumeric);
addParameter(p, 'armijo_c1', 1e-4, @isnumeric);

parse(p, varargin{:});

dobs = p.Results.data;
assert(~isempty(dobs), 'Empty dataset.');

t = p.Results.times;
assert(~isempty(t), 'Empty time range.');

assert(length(t) == length(dobs), 'Data and times not equal.');

r = p.Results.obs;
assert(r > 0.0, 'r must be > 0.');

rho = p.Results.rho;
assert(~isempty(rho), 'Missing starting model layer resistivities.');

thk = p.Results.thickness;
assert(length(thk) == length(rho) - 1, 'Wrong layer thicknesses.');

fix = p.Results.fix;

assert(length(fix) == length(rho), 'Wrong number of fixed params.');
assert(~all(fix), 'At least one parameter must be free.');

if isempty(fix)
    fix = zeros(size(rho));
end

update = ~fix;
P = eye(length(rho));
if ~all(update)
    P(find(~update), :) = []; %#ok
end
WTW = eye(min(size(P')));

verbose = p.Results.verbose;
max_gn_iterations = p.Results.maxit;

lambda = p.Results.lambda;
lambda_min = p.Results.lambda_min;
lambda_scaling_factor = p.Results.reduce_lambda;
lambda_droptol = p.Results.lambda_threshold;
pert = p.Results.pert;
scale_data = p.Results.scale_data;
a = p.Results.scale_asinh;

scalefn = @(x, a) scale_data(x ./ a);

c1 = p.Results.armijo_c1;
step_size_min = p.Results.linesearch_threshold;

tol = p.Results.goal_normgc;
obj_fn_goal = p.Results.goal_obj;
goal_obj_diff = p.Results.goal_obj_diff;

% Some useful function handles
%
getData = @(rho) getVMDLayeredTransient([t(1) t(end)], r, rho, thk, 0.0, 1.0);

residual = @(dobs, d, a) scalefn(dobs, a) - scalefn(d, a);

data_norm = @(r) 0.5 * norm(r)^2;

objective_function = @(r) data_norm(r);

% Transform model parameters, here we take the natural log of rho
%
log_p = log(rho);

% Calculate response of starting model
%
d = getData(rho);

% Calculate data residual for starting model
%
data_residual = residual(dobs, d, a);

obj_fn_current = ...
    data_norm(data_residual);

% Initialize GN monitoring
%
obj_fn = zeros(max_gn_iterations, 1);
lambda_history = zeros(max_gn_iterations, 1);
norm_gc = Inf();
obj_fn_diff = Inf();

gn_it_counter = 0;

obj_fn(1) = obj_fn_current;

% Gauss-Newton iterations start here.
% Iterations will be carried out as long as
% 
% - the norm of the gradient J^T * r is larger than a predefined tolerance
% - maximum iteration is not reached
% - the norm of the objective function is larger than a predefined value
% - the difference between subsequent doesn't fall below a threshold.
%
while (norm_gc > tol) && ...
        (gn_it_counter < max_gn_iterations) && ...
        (obj_fn_current > obj_fn_goal) && ...
        (obj_fn_diff > goal_obj_diff)
    %
    % We have data d, current model log_p, and have evaluated the objective
    % function.
    %
    % Jacobian is defined as J = du/dm * dm / drho.
    %
    J = getJ('data', d, ...
        'rho', exp(log_p), ...
        'thickness', thk, ...
        'times', t, ...
        'obs', r, ...
        'pert', pert, ...
        'scale_asinh', a, ...
        'protem', false);
    
    % Restrict problem to "active" model parameters by only keeping associated columns of J
    %
    J = J * P';
    
    % Gradient J^T * r = -\grad \Phi
    %
    gc = J' * data_residual;
    norm_gc = norm(gc);
    
    % GN search direction, restricted to active parameters
    %
    gn_dir = P' * ((J' * J + lambda * WTW) \ gc);
    
    if verbose
        fprintf('It. %3d: ', gn_it_counter);
        fprintf('JTr: %1.6e ', norm_gc);
    end
    
    % Line Search:
    % Take a fraction of the GN direction and evaluate the response
    %
    step_size = 1.0;
    log_p_proposed = log_p + step_size * gn_dir;
    d_proposed = getData(exp(log_p_proposed));
    
    obj_fn_proposed = objective_function(...
        residual(dobs, d_proposed, a));
    obj_fn_diff = obj_fn_current - obj_fn_proposed;
    
    % Directional derivative for Armijo rule
    % must be negative for descent direction
    dir_deriv = -gn_dir' * (P' * gc);
    
    gn_it_counter = gn_it_counter + 1;
    
    while (obj_fn_diff < -c1 * step_size * dir_deriv)
        step_size = step_size / 2.0;
        if step_size < step_size_min
            error('Line search breakdown. Exit.');
        end
        log_p_proposed = log_p + step_size * gn_dir;
        d_proposed = getData(exp(log_p_proposed));
        obj_fn_proposed = objective_function(...
            residual(dobs, d_proposed, a));
        obj_fn_diff = obj_fn_current - obj_fn_proposed;
    end
    
    obj_rel_drop = obj_fn_diff / obj_fn_current;
    
    % Update parameter after line search
    %
    rho = exp(log_p_proposed);
    log_p = log_p_proposed;
    
    % Update current model response, data residual, and objective function
    %
    d = d_proposed;
    data_residual = residual(dobs, d, a);
    obj_fn_current = obj_fn_proposed;
    
    % Decrease regularization parameter lambda whenever the relative change
    % of the objective function falls below a given drop tolerance.
    % In the first Gauss-Newton iteration, the given value of lambda should
    % always be accepted.
    %
    if obj_rel_drop < lambda_droptol && gn_it_counter > 1
        lambda = lambda * lambda_scaling_factor;
    end
    
    % Reject further decrease of lambda when a given threshold is reached
    %
    if lambda < lambda_min
        lambda = lambda_min;
    end
    
    % Book-keeping
    %
    lambda_history(gn_it_counter) = lambda;
    obj_fn(gn_it_counter) = obj_fn_current;
    
    if verbose
        fprintf('LS: alpha=%1.3e lambda %1.3e obj_fn %1.3e obj_rel_drop %1.3e\n', ...
            step_size, lambda, obj_fn_current, obj_rel_drop);
        
        figure(1);
        loglog(t, abs(dobs), t, abs(d));
        plotTransient(t, dobs, [], 'r');
        hold on
        plotTransient(t, d, [], 'b');
        legend('dobs', '', 'd', '');
        hold off
        set(gcf, 'position', [1 50 560 420]);
        figure(2);
        plot_model(rho, thk);
        set(gcf, 'Position', [613 50 560 419]);
        drawnow();
    end
end

% Return values of model parameters, objective function convergence, and
% change of lambda during iterations
%
m = rho;
lambda_history = lambda_history(1:gn_it_counter);
obj_fn = obj_fn(1:gn_it_counter);

% EOF