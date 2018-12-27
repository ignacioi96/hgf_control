%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants and init values

% length of observed time interval
time_step = 1;
time_interval = 1:time_step:500;


%% Plots
hgf_grapher(@looper);


%% Start of action and belief updates
function [u, mus, x, actions, env_effects, action_effects, S_prediction,...
    S_control, mean_sq_error_control] = looper(time_interval,...
    belief_lambda,belief_alpha,belief_omega,belief_kappa, actual_lambda,...
    actual_alpha,belief_theta, env_effect, mu_des, pi_des, x_init,...
    mu1_init, mu2_init, mu_init_gaussian, pi_init_gaussian,...
    env_effect_func, model_type, action_type, env_effect_period)

    %% initialize all values

    % environment
    % ~~~~~~~~~~
    
    x = x_init * ones(1, time_interval);        % actual value of env 
    u = zeros(1, time_interval);                % perceived value of x
    actions = zeros(1, time_interval);          % agent actions
    action_effects = zeros(1, time_interval);   % effects of agent's actions
    env_effects = zeros(1, time_interval);      % env actions
    
    % internal model of the agent
    % ~~~~~~~~~~

    % hierarchical complexity of the model 
    n_lvls = 2;

    % Estimates of x_1 and x_2 (mu values per level)
    mus = zeros(n_lvls, time_interval);

    % Precision of the estimates of x_1 and x_2 (mu values per level)
    precisions = ones(n_lvls, time_interval);

    % Agents prediction errors on the value of u
    u_pred_errors = zeros(1, time_interval);

    % Volatility pred errors for each level
    volatility_pred_errors = zeros(n_lvls, time_interval);

    if(model_type == "HGF")
        mus(1,:) = mu1_init;
        mus(2,:) = mu2_init;
    elseif(model_type == "simple Gaussian model")
        mus(1,1) = mu_init_gaussian;
        precisions(1) = pi_init_gaussian;
    end
    
    % surprise needs to be initialized
    S_prediction = zeros(1, time_interval);
    S_control = zeros(1, time_interval);
    mean_sq_error_control = zeros(1, time_interval);
    
    %% belief update and actions taken
    % This is the core structure of the program.
    for i=2:time_interval+1
           
        % generate a new sensation
        u(i) = sampleU(x(i-1), actual_alpha);

        if(model_type == "HGF")
            % update the internal model
            [muhat, pihat, u_pred_errors(i),...
                mus(:,i), precisions(:,i),...
                volatility_pred_errors(:,i)] = hgf(u(i), mus(:,i-1),...
                precisions(:,i-1), actions(i-1), belief_lambda,...
                belief_alpha, belief_omega, belief_kappa, belief_theta);
        elseif(model_type == "simple Gaussian model")
            [muhat, pihat, u_pred_errors(i),...
                mus(:,i), precisions(:,i),...
                volatility_pred_errors(:,i)] = gaussian(u(i),...
                mus(:, i-1), precisions(:, i-1), actions(i-1),...
                belief_lambda, belief_alpha); 
        end

        % calculate actions
        % Notice that this is not using the same prediction error as the
        % hgf function. The model has already been updated here.
        if(action_type == "model-based control")
            actions(i) = act_model(mu_des, pi_des, mus(1,i),...
                precisions(1,i));
        elseif(action_type == "homeostatic control")
            actions(i) = act_homeostasis(mu_des, pi_des, u(i));
        end
        
        % use actions and external influences to change the environment
        [x(i), action_effects(i-1)] = changeEnv_action(...
            actions(i), x(i-1), actual_lambda);
        
        % calculate predictive surprise and control surprise
        S_prediction(i) = -0.5*(log(pihat(1))-pihat(1)*u_pred_errors(i)^2);
        S_control(i) = -0.5*(log(pi_des(1))-pi_des(1)*(u(i)-mu_des)^2);
        mean_sq_error_control(i) = (x(i)-mu_des)^2;
        
        [x(i), env_effects(i-1)] = changeEnv_envEff(i, x(i),...
            env_effect, env_effect_func, env_effect_period);
     
    end
end

%% Learning

% Gaussian models
function [muhat, pihat, dau,...
    mu, precision, da] = gaussian(u, mu, precision, action, lambda, alpha)
    
    % input format:
    % mu, precision, u, action are scalars (we sample stepwise)
    
    muhat = NaN(2,1);
    pihat = NaN(2,1);
    
    muhat(1) = mu(1) + lambda*action;
    pihat(1) = precision(1);
    
    % Input/Value prediction error
    dau = u-g(muhat(1));
    
    precision(1) = pihat(1) + 1/alpha;    
    mu(1) = muhat(1) + ((1/alpha)/(precision(1)))*(u-g(muhat(1)));
    
    da = NaN(2,1);
end

% HGF
function [muhat, pihat, dau,...
    mu, precision, da] = hgf(u, mu, precision, action, lambda, alpha,...
    omega, kappa, theta)
    
    % input format:
    % mu, precision are 2x1
    % u, action are scalars (we sample stepwise)
    
    % reading out the parameters
    th = exp(theta);
    al = alpha;
    ka = kappa;
    om = omega;
    
    % Lilian told us we don't need these:
    t = 1;
    rho = 0;
    
    % initialization
    da = zeros(2,1);
    muhat = zeros(length(mu),1);
    pihat = zeros(length(precision),1);
    

    % First Level
    % ~~~~~~~~~~~
    
    % Prediction
    % -> action is included here now
    muhat(1) = mu(1) + t*rho + lambda*action;

    % Precision of prediction
    pihat(1) = 1/(1/precision(1) + t*exp(ka*mu(2)+om));

    % Input/Value prediction error
    dau = u-muhat(1);

    % Updates
    precision(1) = pihat(1) + 1/al;
    mu(1) = muhat(1) + 1/pihat(1) * 1/(1/pihat(1) + al) * dau;

    % Volatility prediction error
    da(1) = (1/precision(1) + (mu(1)-muhat(1))^2) *pihat(1)-1;
    
    
    % Last (=second) level: volatility
    % Here we might encounter problems because we also learn the volatility
    % that is induced by the agents actions.
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
    mu(2) = muhat(2) +1/2 *1/precision(2) *ka(1) *w(1) *da(1);

    % Volatility prediction error
    da(2) = (1/precision(2) +(mu(2)-muhat(2))^2) *pihat(2) -1;
end

%% Action
% calculate action based on internal model and desired state
function a = act_model(mu_des, pi_des, mu_1, pi_1)
    a = pi_des*(mu_des - mu_1);
    % = precision of homeostatic belief * prediction error of model
end

% calculate action solely based on input
function a = act_homeostasis(mu_des, pi_des, u)
    a = pi_des*(mu_des - u);
    % = precision of homeostatic belief * prediction error
end

% effector function: scaling by an efficacy factor
function effect = f(action, lambda)
    effect = lambda*action;
end

%% updating x
% action effects taken into account
function [x_new, action_effect] = changeEnv_action(action, x, lambda)
    
    % return the external factor
    action_effect = f(action, lambda);
    
    % update
    x_new = x + action_effect;
end

% environment effects taken into account
function [x_new, env_effect] = changeEnv_envEff(time_point,...
    x, env_effect, func, env_effect_dilation)
    
    % calculate the external perturbation
    % func basically solves as the derivative of the real effect here
    external_factor = env_effect*func(time_point/env_effect_dilation);
    
    % return the external factor
    env_effect = external_factor;
    
    % update
    x_new = x + env_effect;
end

%% "generates sensations"
% samples u from the probabilistic sensor, given the true value of x
function u = sampleU(actual_mean, var_sensor)
    u = normrnd(actual_mean, sqrt(var_sensor));
    % write g(actual_mean) to include the sensor function
end

%% sensor function
% This part became irrelevant when using the HGF, because there an identity
% sensor function (+noise of course) is implicitly assumed.
function u = g(x)
    u = x;
end

% derivative of sensor function
function u = dg(x)
    u = 1;
end
