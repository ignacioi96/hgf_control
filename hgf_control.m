% Free parameters that can be changed:
% ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set constants

% data does not temporarily evolve
x = 5;

% precision of the sensor
sigma_data = 100;

% belief of perceived reality is a normal distribution
mu_per = 0;
sigma_per = 0.0;

% desired reality 
mu_des = 0;
sigma_des = 0.01;

% speed of actions on environment
lambda = 1;
lambda_inv = 1/lambda;

time_step = 1;
time_interval = 1:time_step:100;

x_val = time_interval;
action_timeline = time_interval;

density_interval = -10:.05:10;
x_perceived = zeros(length(time_interval), length(density_interval));

per_mus = zeros(length(time_interval),1);
per_sigmas = zeros(length(time_interval),1);

val_pred_errors = zeros(length(x_perceived),1);
volatility_pred_errors = zeros(length(x_perceived),1);

n_lvls = 2;
mus = zeros(n_lvls,1);
precisions = zeros(n_lvls,1);

for i=time_interval
    x_perceived(i,:) = normpdf(density_interval,mu_per,sqrt(1/sigma_per));

    y_t = sampleY(x, sigma_data);
%     [mu_per, sigma_per] = update(mu_per, sigma_per, y_t, sigma_data);
 
    %action
    x_val(time_interval==i) = x;

    per_mus(i) = mu_per;
    per_sigmas(i) = sigma_h;
    
    % hgf(vals, mu, per)
    [mu_per, sigma_per, dau, mu, sigma, da] = hgf(x_val,...
        mu_per, sigma_per);
    val_pred_errors(i);
    
    action_timeline(time_interval==i) = action(mu_des,...
       sigma_des, sampleY(x, sigma_data), x); 
   
    x = changeEnv(time_step, lambda, mu_des,...
        sigma_des, y_t, x);%sampleY(x, sigma_data), x);
    
%     tapas_fitModel(mu_per,x_val);

end
% response mosel unitsq 
p1 = subplot(1,2,1);
plot(x_val);
hold on;

colormap(p1, winter);
plot(per_mus(:,1));
hold on;

mean_val = mean(x_val);
plot(mean_val*ones(size(time_interval)));
axis square;
title('Real X');
legend('Real', 'Perceived', 'Mean of Real Value');

p3 = subplot(1,2,2);
plot(time_interval, action_timeline*0.2);
colormap(p3,spring);
axis square;
title('Actions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hgf
function [muhat, pihat, dau, mu, pre, da] = hgf(val, mu, pre)
   % inputs are 2 d, 1 input per level

    % First Level
    % ~~~~~~~~~~~
    l_val = length(val);
    % no idea what this are need to ask
    t = 0.5;
    rho = 1;
    ka = 0.5;
    om = .01;
    th = exp(0.5);
    al = 1;
    da = 0;
    
    % Prediction
    muhat(1) = mu(1) + t*rho;

    % Precision of prediction
    pihat(1) = 1/(1/pre(1)) + t *exp(ka*mu(1)+om);

    % Input/Value prediction error
    dau = val-muhat(1);

    % Updates
    pre = pihat(1) + 1/al;
    mu = muhat(1) + 1/pihat(1) * 1/(1/pihat(1) + al) * dau;

    % Volatility prediction error
    da = (1/pre(1) + mu(1)-muhat(1)^2)*pihat(1)-1;
    
    
    % Last level (we only have 2 levels later this will be volatility)
    % right now this is action influenced variable observation
    % ~~~~~~~~~~
    % Prediction
    muhat(2) = mu(2) +t*rho;

    % Precision of prediction
    pihat(2) = 1/(1/pre(2) +t*th);

    % Weighting factor
    v(2)   = t *th;
    v(1) = t *exp(ka(1) *mu(2) +om(1));
    w(1) = v(1) *pihat(1);

    % Updates
    pre(2) = pihat(2) +1/2 *ka(1)^2 *w(1) *(w(1) +(2 *w(1) -1) *da(1));

%     if pi(k,l) <= 0
%         error('tapas:hgf:NegPostPrec', 'Negative posterior precision. Parameters are in a region where model assumptions are violated.');
%     end

    mu(2) = muhat(2) +1/2 *1/pre(2) *ka(1) *w(1) *da(1);

    % Volatility prediction error
    da(2) = (1/pre(2) +(mu(2)-muhat(2))^2) *pihat(2) -1;

end

% effect of action on environment
function x = changeEnv(time_step, lambda, mu_prior, sigma_prior, y_t, x)
    x = x + time_step*(1/lambda)*f(action(mu_prior, sigma_prior, y_t, x));
end

% effector function
function effect = f(action)
    effect = action;
end

% what does the agent do?
function a = action(mu_prior, sigma_prior, y_t, x)
    a = -sigma_prior*(y_t-g(mu_prior))*dg(x);
end
% updating beliefs based on data
% hgf outputs prediction errors
% response models
% noise param on action simulation (otherwise need to be super precise
% about )

% generates sensations
function y = sampleY(mean, sigma_data)
    y = normrnd(g(mean), sqrt(1/sigma_data));
end

% sensor function
function y_mean = g(x)
    y_mean = x;
end

% derivative of sensor function
function value = dg(x)
    value = 1;
end