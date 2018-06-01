% Free parameters that can be changed:
% ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and init values

% length of observed time interval
time_step = 1;
time_interval = 1:time_step:500;

% variance of the sensor of the agent
alpha = 0.05;

%% Plots
hgf_grapher(@looper);

%% Start of action and belief updates
function [u, mus, x, actions] = looper(time_interval, belief_lambda,...
    belief_alpha,belief_omega,belief_kappa, actual_lambda, actual_alpha,...
    belief_theta, env_effect, mu_des, pi_des, x_init)

    %% X, values of environment
    % Actual value of env 
    x = zeros(1, time_interval);
    x(:) = x_init;
    % Perceived value of x
    u = zeros(1,time_interval);
    
    %% Belief model that the agent has of the env
    % hierarchical complexity of the model 
    n_lvls = 2;

    % Estimates of x_1 and x_2 (mu values per level)
    mus = ones(n_lvls,time_interval);

    % Precision of the estimates of x_1 and x_2 (mu values per level)
    precisions = ones(n_lvls,time_interval);

    % Agents prediction errors on the value of u
    u_pred_errors = zeros(1,time_interval);

    % Volatility pred errors for each level
    volatility_pred_errors = zeros(n_lvls,time_interval);

%     %% Desired model of env
%     mu_des = 0;
%     pi_des = 0.01;

    %% Actions performed by agent on x 
    actions = zeros(1,time_interval);
    
    %% belief update and actions taken
    for i=2:time_interval
        u(i) = sampleU(x(i-1), actual_alpha);

        [muhat, prehat, u_pred_errors(i),...
            mus(:,i), precisions(:,i),...
            volatility_pred_errors(:,i)] = hgf(u(i), mus(:,i-1),...
            precisions(:,i-1), actions(i-1), belief_lambda,...
            belief_alpha, belief_omega, belief_kappa, belief_theta);

        actions(i) = act(mu_des, pi_des, mus(1,i),...
            precisions(1,i), actual_lambda); 

        x(i) = changeEnv(i, actions(i), x(i-1), env_effect, actual_lambda);
    end
end

%% HGF
function [muhat, pihat, dau,...
    mu, precision, da] = hgf(u, mu, precision, action, lambda, alpha,...
    omega, kappa, theta)
   % inputs are 2 d, 1 input per level

    % First Level
    % ~~~~~~~~~~~
    length_u = length(u);
    % no idea what this are need to ask
    t = 1;
    rho = 0;
    th = exp(theta);
    al = alpha;
    ka = kappa;
    om = omega;
    da = 0;
    muhat = zeros(length(mu),1);
    pihat = zeros(length(precision),1);
    

    % Prediction
    muhat(1) = mu(1) + t*rho + lambda*action;

    % Precision of prediction
    pihat(1) = 1/(1/precision(1) + t *exp(ka*mu(2)+om));

    % Input/Value prediction error
    dau = u-muhat(1);

    % Updates
    precision(1) = pihat(1) + 1/al;
    mu(1) = muhat(1) + 1/pihat(1) * 1/(1/pihat(1) + al) * dau;

    % Volatility prediction error
    da(1) = (1/precision(1) + (mu(1)-muhat(1))^2) *pihat(1)-1;
    
    
    % Last level (we only have 2 levels later this will be volatility)
    % right now this is action influenced variable observation
    % ~~~~~~~~~~
    % Prediction
    muhat(2) = mu(2) +t*rho;

    % Precision of prediction
    pihat(2) = 1/(1/precision(2) +t*th);

    % Weighting factor
    v(2)   = t *th;
    v(1) = t * exp(ka(1) *mu(2) +om(1));
    w(1) = v(1) * pihat(1);

    % Updates
    precision(2) = pihat(2) +1/2 *ka(1)^2 *w(1) *(w(1) +(2 *w(1) -1) *da(1));

%     if pi(k,l) <= 0
%         error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
%     end

    mu(2) = muhat(2) +1/2 *1/precision(2) *ka(1) *w(1) *da(1);

    % Volatility prediction error
    da(2) = (1/precision(2) +(mu(2)-muhat(2))^2) *pihat(2) -1;
end

% effect of action on environment
function x_new = changeEnv(time_point, action, x, env_effect, lambda)
    external_factor = env_effect*sin(0.01*time_point);
    x_new = x + f(action, lambda) + external_factor;
end

% effector function
function effect = f(action, lambda)
    effect = lambda*action;
end

%% Action
% based on 
% what does the agent do?
function a = act(mu_des, pi_des, mu_1, pi_1, lambda)
    a = mu_des - (1/sqrt(pi_1))*mu_1;
end
% updating beliefs based on data
% hgf outputs prediction errors
% response models
% noise param on action simulation (otherwise need to be super precise
% about )

% generates sensations
function y = sampleU(mean, var_data)
    y = normrnd(g(mean), sqrt(var_data));
end

% sensor function
function y_mean = g(x)
    y_mean = x;
end

% derivative of sensor function
function value = dg(x)
    value = 1;
end