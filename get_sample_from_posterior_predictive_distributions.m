function [V_posterior_predictive_1, ...
    V_posterior_predictive_2, ...
    V_posterior_predictive_sum] = ... 
    get_sample_from_posterior_predictive_distributions(experiment_to_consider, ...
    rc_sampled_from_posterior,...
Kc_sampled_from_posterior, ...
lc_sampled_from_posterior,...
rr_sampled_from_posterior,...
Kr_sampled_from_posterior,...
lr_sampled_from_posterior,...
V0_sum_sampled_from_posterior, sigma_sampled_from_posterior)
% Gets a sample of parameters from posterior, (rr, Kr, lr, rc, Kc, lc, V0, sigma)
% and returns posterior predictive distributions using these parameters,
% solving the ODE and returning sampled Vc and Vr values drawn from a
% normal distribution centred at each time point around the Vc or Vr
% solution of the ODE
%
% For the input parameter experiment_to_consider: its value is between 1 and 5, 
% and each value encodes and maps to a particular ratio chosen in a certain experiment.
% The mapping is the following :
% 1 -> 1:0
% 2 -> 3:1
% 3 -> 1:1
% 4 -> 1:3
% 5 -> 0:1

%% 0. Parameters setting

% Enter parameters

odds_of_initial_values = [1,0;
                          3,1;
                          1,1;
                          1,3;
                          0,1];

beginning_time = 1;
end_time = 15;

time_first_individual_data = 10;
sampling_timestep_individual_data = 5;
number_of_experimental_repeats = 12;


% Calculations from parameters

    % Sampled times
offset_time_first_individual_data = time_first_individual_data - beginning_time;
sampled_times_individual_data = beginning_time+offset_time_first_individual_data:...
    sampling_timestep_individual_data:...
    end_time;
sampled_times_sum_data = [3, 5, 7, 10, 13, 15];

    % Ratio of odds to get V0 of individual dependent variables (Vc and Vr)
sum_of_odds = odds_of_initial_values(experiment_to_consider,1) + ...
odds_of_initial_values(experiment_to_consider,2);
ratio_V0_C = odds_of_initial_values(experiment_to_consider,1) / ...
sum_of_odds;
ratio_V0_R = odds_of_initial_values(experiment_to_consider,2) / ...
sum_of_odds;

V0_sampled_from_posterior = [V0_sum_sampled_from_posterior*ratio_V0_C, ...
    V0_sum_sampled_from_posterior*ratio_V0_R]; 

    % Miscellaneous
tspan = [beginning_time , end_time];

sigma_sum_normal_sampled_from_posterior = sigma_sampled_from_posterior / ...
    sqrt(number_of_experimental_repeats);

%% 1. solve the ODE model with the previously sampled (rr, Kr, lr, rc, Kc, lc, V0) 
% and a defined ratio r (1:0, 1:1, 3:1...), the result are mu_R, mu_C, mu_sum

[t,V_posterior_predictive] = ode45(...
    @(t,V) ode_model(t, V, ...
    rc_sampled_from_posterior, ...
    Kc_sampled_from_posterior, ...
    lc_sampled_from_posterior, ...
    rr_sampled_from_posterior, ...
    Kr_sampled_from_posterior, ...
    lr_sampled_from_posterior), ...
        tspan,...
        V0_sampled_from_posterior);

    
%% 2. draw samples from Normal(mu_R (== V_posterior_predictive(:,1)), sigma),
% Normal(mu_C (== V_posterior_predictive(:,2)), sigma), 
% Normal(mu_sum, sigma).

V_posterior_predictive_1 = normrnd(V_posterior_predictive(:,1),sigma_sampled_from_posterior);
V_posterior_predictive_2 = normrnd(V_posterior_predictive(:,2),sigma_sampled_from_posterior);
V_posterior_predictive_sum = normrnd(V_posterior_predictive(:,1) + ...
    V_posterior_predictive(:,2),sigma_sum_normal_sampled_from_posterior);

V_posterior_predictive_1 = interp1(t,V_posterior_predictive_1,sampled_times_individual_data);
V_posterior_predictive_2 = interp1(t,V_posterior_predictive_2,sampled_times_individual_data);
V_posterior_predictive_sum = interp1(t,V_posterior_predictive_sum, sampled_times_sum_data);




end

