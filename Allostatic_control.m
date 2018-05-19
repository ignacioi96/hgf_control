% Free parameters that can be changed:
% ?



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set constants

% data does not temporarily evolve
x = 5;

% noise in the sensor
pi_data = 0.5;

% belief is a normal distribution
mu_t = 0;
pi_t = 1;

% speed of actions on environment
lambda = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST of action

time_step = 0.01;
grid = 0:time_step:8;
model = grid;
action_timeline = grid;
for i=grid
    model(grid==i) = x;
    action_timeline(grid==i) = action(mu_t,...
        pi_t, sampleY(x, pi_data), x);
    
    x = changeEnv(time_step, lambda, mu_t,...
        pi_t, sampleY(x, pi_data), x);
end
plot(grid, model)
hold on
plot(grid, action_timeline.*0.3)




% TEST of belief update

% grid = -3:.01:8;
% model = normpdf(grid,mu_t,sqrt(1/pi_t));
% plot(grid, model)
% 
% for i=1:50
%     hold on
%     model = normpdf(grid,mu_t,sqrt(1/pi_t));
%     plot(grid, model)
% 
%     y_t = sampleY(x, pi_data);
%     [mu_t, pi_t] = update(mu_t, pi_t, y_t, pi_data);
% end




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