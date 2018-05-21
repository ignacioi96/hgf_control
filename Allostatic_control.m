% Free parameters that can be changed:
% ?



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set constants

% data does not temporarily evolve
x = 5;

% noise in the sensor
pi_data = 0.05;

% belief of perceived reality is a normal distribution
mu_per = 0;
pi_per = 1;

% desired reality 
mu_des = 5;
pi_des = 2;

% speed of actions on environment
lambda = 1;
lambda_inv = 1/lambda;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST of action

% time_step = 1;
% time_interval = 1:time_step:10;
% 
% x_val = time_interval;
% action_timeline = time_interval;
% for i=time_interval
%     x_val(time_interval==i) = x;
%     action_timeline(time_interval==i) = action(mu_des,...
%         pi_des, sampleY(x, pi_data), x);
%     
%     x = changeEnv(time_step, lambda, mu_des,...
%         pi_des, sampleY(x, pi_data), x);
% end
% plot(time_interval, x_val)
% hold on
% plot(time_interval, action_timeline.*0.2)




% TEST of perceived reality update

% density_interval = -3:.01:8;
% x_per = zeros(length(time_interval), length(density_interval));
% 
% for i=time_interval
%     x_per(i,:) = normpdf(density_interval,mu_per,sqrt(1/pi_per));
% 
%     y_t = sampleY(x, pi_data);
%     [mu_per, pi_per] = update(mu_per, pi_per, y_t, pi_data);
% end
% surf(density_interval.',time_interval, x_per)



% TEST both

% for i=time_interval

% end
% plot(time_interval, x_val)
% hold on
% plot(time_interval, action_timeline.*0.2)

time_step = 1;
time_interval = 1:time_step:100;

x_val = time_interval;
action_timeline = time_interval;

density_interval = -3:.01:8;
x_per = zeros(length(time_interval), length(density_interval));

for i=time_interval
    x_per(i,:) = normpdf(density_interval,mu_per,sqrt(1/pi_per));

    y_t = sampleY(x, pi_data);
    [mu_per, pi_per] = update(mu_per, pi_per, y_t, pi_data);
    
% %     action
    x_val(time_interval==i) = x;
    action_timeline(time_interval==i) = action(mu_des,...
        pi_des, sampleY(x, pi_data), x);
    
    x = changeEnv(time_step, lambda, mu_des,...
        pi_des, sampleY(x, pi_data), x);

end
surf(density_interval.',time_interval, x_per)




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
% TODO: derive the update equation here
% TODO: combine this with control
% Is this the same if we pass multiple data points at once?
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