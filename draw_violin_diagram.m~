function draw_violin_diagram()

experiment_to_consider = 1;
number_of_posterior_predictive_samples = 10;

%% Data

V_sum_data_average_PC3_0Gy = [0.0636679167 0.0661274167 0.0987616667 0.16010675 0.25829375 0.280966;
    0.0574364167 0.0799836667 0.1146740833 0.2577545833 0.41375175 0.4119126667;
    0.0589646667 0.0776658333 0.1201630833 0.27822025 0.5145295 0.557941;
    0.0544568333 0.08331 0.1326801667 0.35862575 0.5359255 0.523286;
    0.0503110833 0.0921101667 0.1695063333 0.4100345833	0.5936558333 0.6179516667];

V_individual_control_data_PC3_0Gy = [0.99 0.98;
    0.56 0.45;
    0.32 0.18;
    0.12 0.06;
    0.00 0.00];

V_individual_resistant_data_PC3_0Gy = [0.00 0.00;
    0.43 0.53;
    0.67 0.80;
    0.87 0.92;
    0.99 0.97];

%% Ugly recalculation of number_of_sampled_times_individual_variables 
% and number_of_sampled_times_sum

beginning_time = 1;
end_time = 15;

time_first_individual_data = 10;
sampling_timestep_individual_data = 5;

offset_time_first_individual_data = time_first_individual_data - beginning_time;
sampled_times_individual_data = beginning_time+offset_time_first_individual_data:...
    sampling_timestep_individual_data:...
    end_time;
sampled_times_sum_data = [3, 5, 7, 10, 13, 15];

number_of_sampled_times_individual_variables = size(sampled_times_individual_data,2);
number_of_sampled_times_sum = size(sampled_times_sum_data,2);

%% The code really starts here

V_posterior_predictive_1_array = nan(number_of_posterior_predictive_samples, ...
    number_of_sampled_times_individual_variables);
V_posterior_predictive_2_array = nan(number_of_posterior_predictive_samples, ...
    number_of_sampled_times_individual_variables);
V_posterior_predictive_sum_array = nan(number_of_posterior_predictive_samples, ...
    number_of_sampled_times_sum);

for posterior_predictive_sample_number = 1 : number_of_posterior_predictive_samples
    %% Choose method of getting posterior parameter values
    [rc_sampled_from_posterior,...
    Kc_sampled_from_posterior, ...
    lc_sampled_from_posterior,...
    rr_sampled_from_posterior,...
    Kr_sampled_from_posterior,...
    lr_sampled_from_posterior,...
    V0_sum_sampled_from_posterior, sigma_sampled_from_posterior] = ABC_normal_prior_fitted_sigma();
    
    %% Plug sampled parameter values in get_sample_from_posterior_predictive_distributions
    [V_posterior_predictive_1, ...
        V_posterior_predictive_2, ...
        V_posterior_predictive_sum] = ... 
        get_sample_from_posterior_predictive_distributions(experiment_to_consider, ...
        rc_sampled_from_posterior,...
    Kc_sampled_from_posterior, ...figure;
distributionPlot(V_posterior_predictive_sum_array,'xValues',sampled_times_sum_data);
plotSpread(V_sum_data_average_PC3_0Gy(experiment_to_consider,:),'xValues',sampled_times_sum_data);
    lc_sampled_from_posterior,...
    rr_sampled_from_posterior,...
    Kr_sampled_from_posterior,...
    lr_sampled_from_posterior,...
    V0_sum_sampled_from_posterior, sigma_sampled_from_posterior);

    %% Store sampled posterior predictives in 
    
    V_posterior_predictive_1_array(posterior_predictive_sample_number, :) = ...
        V_posterior_predictive_1;
    
    V_posterior_predictive_2_array(posterior_predictive_sample_number, :) = ...
        V_posterior_predictive_2;
    
    V_posterior_predictive_sum_array(posterior_predictive_sample_number, :) = ...
        V_posterior_predictive_sum;
end

figure;
distributionPlot(V_posterior_predictive_sum_array,'xValues',sampled_times_sum_data);
plotSpread(V_sum_data_average_PC3_0Gy(experiment_to_consider,:),'xValues',sampled_times_sum_data);


end

