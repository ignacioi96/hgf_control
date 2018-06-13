% Free parameters that can be changed:
% ?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set constants

% data does not temporarily evolve
x = 5;

% precision of the sensor
pi_data = 0.1;

% belief of perceived reality is a normal distribution
mu_per = 0;
pi_per = 0.0;

% desired reality 
mu_des = 0;
pi_des = 0.01;

% speed of actions on environment
lambda = 1;
lambda_inv = 1/lambda;

% initializing arrays to store values in loop
time_step = 1;
time_interval = 1:time_step:300;

x_val = time_interval;
action_timeline = time_interval;

density_interval = -10:.05:10;
x_per = zeros(length(time_interval), length(density_interval));
x_val_mat = x_per;
mus = zeros(1,length(time_interval));

for i=time_interval
    x_per(i,:) = normpdf(density_interval,mu_per,sqrt(1/pi_per));

    y_t = sampleY(x, pi_data);
    [mu_per, pi_per] = update(mu_per, pi_per, y_t, pi_data);
    mus(i) = mu_per;
    
    %action
    x_val(time_interval==i) = x;
    x_val_mat(i, round((x+10)*20)) = 1;
    action_timeline(time_interval==i) = action(mu_des,...
       pi_des, sampleY(x, pi_data), x);
    
    x = changeEnv(time_step, lambda, mu_des,...
        pi_des, y_t, x);%sampleY(x, pi_data), x);

end

p1 = subplot(1,2,1);
plot(x_val);
hold on;
colormap(p1, winter);
plot(mus);
hold on;
mean_val = mean(x_val);
plot(mean_val*ones(size(time_interval)));
axis square;
title('Values of X');
legend('Real', 'Believed', 'Mean of Real Value');

p3 = subplot(1,2,2);
plot(time_interval, action_timeline);
colormap(p3,spring);
axis square;
title('Actions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% effect of action on environment
function x = changeEnv(time_step, lambda, mu_prior, pi_prior, y_t, x)
    x = x + time_step*(1/lambda)*f(action(mu_prior, pi_prior, y_t, x));
end

% effector function
function effect = f(action)
    effect = action;
end

% what does the agent do?
function a = action(mu_prior, pi_prior, y_t, x)
    a = -pi_prior*(y_t-g(mu_prior))*dg(x);
end
% updating beliefs based on data
% why can we pass pi_data? How do we know that?

function [mu_tplus1, pi_tplus1] = update(mu_t, pi_t, y_t, pi_data)
    pi_tplus1 = pi_t + pi_data;    
    mu_tplus1 = mu_t + (pi_data/(pi_tplus1))*(y_t-g(mu_t));
end

% generates sensations
function y = sampleY(mean, pi_data)
    y = normrnd(g(mean), sqrt(1/pi_data));
end

% sensor function
function y_mean = g(x)
    y_mean = x;
end

% derivative of sensor function
function value = dg(x)
    value = 1;
end