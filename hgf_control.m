% Free parameters that can be changed:
% ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and init values

% length of observed time interval
time_step = 1;
time_interval = 1:time_step:100;

% variance of the sensor of the agent
alpha = 0.8;

% strength of action on the env
lambda = 0.1;
lambda_inv = 1/lambda;

%% X, values of environment

% Actual value of env 
x = zeros(1, length(time_interval));
x(1) = 5;

% Perceived value of x
u = zeros(1,length(time_interval));

%% Belief model that the agent has of the env
% hierarchical complexity of the model 
n_lvls = 2;

% Estimates of x_1 and x_2 (mu values per level)
mus = ones(n_lvls,length(time_interval));

% Precision of the estimates of x_1 and x_2 (mu values per level)
precisions = ones(n_lvls,length(time_interval));

% Agents prediction errors on the value of u
u_pred_errors = zeros(1,length(time_interval));

% Volatility pred errors for each level
volatility_pred_errors = zeros(n_lvls,length(time_interval));

%% Desired model of env
mu_des = 0;
pi_des = 0.01;

%% Actions performed by agent on x 
actions = zeros(1,length(time_interval));

%% Start of action and belief updates
for i=2:length(time_interval)
    u(i) = sampleU(x(i-1), alpha);

    [muhat, prehat, u_pred_errors(i),...
        mus(:,i), precisions(:,i),...
        volatility_pred_errors(:,i)] = hgf(u(i), mus(:,i-1), precisions(:,i-1));
    
    actions(i) = act(mu_des, pi_des, u(i), x(i-1)); 
   
    x(i) = changeEnv(time_step, lambda, actions(i), x(i-1));
end

%% Plots
p1 = subplot(1,2,1);
% Real Value of X
plot(x);
hold on;
colormap(p1, winter);
% Perceived value of X
plot(u);
hold on;
% Mean value of X
mean_val = mean(x);
plot(mean_val*ones(size(time_interval)));
axis square;
title('Real X');
legend('Real', 'Perceived', 'Mean of Real Value');

p2 = subplot(1,2,2);
plot(time_interval, actions);
colormap(p2,spring);
axis square;
title('Actions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% hgf
function [muhat, pihat, dau, mu, precision, da] = hgf(val, mu, precision)
   % inputs are 2 d, 1 input per level

    % First Level
    % ~~~~~~~~~~~
    l_val = length(val);
    % no idea what this are need to ask
    t = 1;
    rho = 0;
    ka = 0.5;
    om = .01;
    th = exp(0.5);
    al = 1;
    da = 0;

    muhat = zeros(length(mu),1);
    pihat = zeros(length(precision),1);
    

    % Prediction
    muhat(1) = mu(1) + t*rho;

    % Precision of prediction
    pihat(1) = 1/(1/precision(1)) + t *exp(ka*mu(1)+om);

    % Input/Value prediction error
    dau = val-muhat(1);

    % Updates
    precision(1) = pihat(1) + 1/al;
    mu(1) = muhat(1) + 1/pihat(1) * 1/(1/pihat(1) + al) * dau;

    % Volatility prediction error
    da(1) = (1/precision(1) + mu(1)-muhat(1)^2)*pihat(1)-1;
    
    
    % Last level (we only have 2 levels later this will be volatility)
    % right now this is action influenced variable observation
    % ~~~~~~~~~~
    % Prediction
    muhat(2) = mu(2) +t*rho;

    % Precision of prediction
    pihat(2) = 1/(1/precision(2) +t*th);

    % Weighting factor
    v(2)   = t *th;
    v(1) = t *exp(ka(1) *mu(2) +om(1));
    w(1) = v(1) *pihat(1);

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
function x_new = changeEnv(time_step, lambda, action, x)
    x_new = x + time_step*(1/lambda)*f(action);
end

% effector function
function effect = f(action)
    effect = action;
end

% what does the agent do?
function a = act(mu_prior, sigma_prior, y_t, x)
    a = -sigma_prior*(y_t-g(mu_prior))*dg(x);
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