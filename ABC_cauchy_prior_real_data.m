function ABC()

number_of_iterations = 3000;
beginning_time = 1;
end_time = 15;
time_first_individual_data = 10;
time_first_sum_data = 3;
sampling_timestep_individual_data = 5;
sampling_timestep_sum_data = 2;
number_of_parameters_to_estimate = 7;
number_of_individual_variables = 2;
number_of_different_ratios_tried = 5;
odds_of_initial_values = [1,0;
                          3,1;
                          1,1;
                          1,3;
                          0,1];
                      

offset_time_first_individual_data = time_first_individual_data - beginning_time;
offset_time_first_sum_data = time_first_sum_data - beginning_time;
tspan = [beginning_time , end_time];

sampled_times_individual_data = beginning_time+offset_time_first_individual_data:...
    sampling_timestep_individual_data:...
    end_time;
sampled_times_sum_data = beginning_time+offset_time_first_sum_data:...
    sampling_timestep_sum_data:end_time;


number_of_sampled_times_individual_variables = size(sampled_times_individual_data,2);
number_of_sampled_times_sum = size(sampled_times_sum_data,2);

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

rc_true_PC3_0Gy=0.236;
Kc_true_PC3_0Gy=0.473;
lc_true_PC3_0Gy=0.655;
rr_true_PC3_0Gy=0.401;
Kr_true_PC3_0Gy=0.773;
lr_true_PC3_0Gy=0.480;
V0_sum_true_PC3_0Gy = 0.031 +  0.015;

total_number_of_sampled_times = (number_of_individual_variables *...
    number_of_sampled_times_individual_variables + ...
    number_of_sampled_times_sum) * ...
    number_of_different_ratios_tried;

tolerance = 0.1*total_number_of_sampled_times;
standard_deviation_noise = 0.01;
number_of_experimental_repeats = 12;
standard_deviation_sum_normal = standard_deviation_noise / sqrt(number_of_experimental_repeats);

accepted_rc_array = nan(number_of_iterations,1);
accepted_Kc_array = nan(number_of_iterations,1);
accepted_lc_array = nan(number_of_iterations,1);
accepted_rr_array = nan(number_of_iterations,1);
accepted_Kr_array = nan(number_of_iterations,1);
accepted_lr_array = nan(number_of_iterations,1);
accepted_V0_sum_array = nan(number_of_iterations,1);

% Generating Cauchy prior
% Initial values for PC3	0Gy
% cauchy_prior_rc = makedist('tLocationScale','mu',0.236,'sigma',1,'nu',1);
% cauchy_prior_Kc = makedist('tLocationScale','mu',0.473,'sigma',1,'nu',1);
% cauchy_prior_lc = makedist('tLocationScale','mu',0.655,'sigma',1,'nu',1);
% cauchy_prior_rr = makedist('tLocationScale','mu',0.401,'sigma',1,'nu',1);
% cauchy_prior_Kr = makedist('tLocationScale','mu',0.773,'sigma',1,'nu',1);
% cauchy_prior_lr = makedist('tLocationScale','mu',0.480,'sigma',1,'nu',1);
% cauchy_prior_V0_sum = makedist('tLocationScale','mu',0.031 + 0.015,'sigma',1,'nu',1);

% another prior, more realistic
cauchy_prior_rc = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_Kc = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_lc = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_rr = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_Kr = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_lr = makedist('tLocationScale','mu',0.5,'sigma',0.5,'nu',1);
cauchy_prior_V0_sum = makedist('tLocationScale','mu',0.05,'sigma',0.5,'nu',1);

% % Parameters for normal pr3,1;
% mean_rc = 2;
% standard_deviation_rc = 1;
% 
% mean_Kc = 2;
% standard_deviation_Kc = 1;
% 
% mean_lc = 2;
% standard_deviation_lc = 1;(iteration_with_certain_ratio)
% 
% mean_rr = 2;
% standard_deviation_rr = 1;
% 
% mean_Kr = 2;
% standard_deviation_Kr = 1;
% 
% mean_lr = 2;
% standard_deviation_lr = 1;

% [V_fakedata_sum,V_fakedata_1, V_fakedata_2,rc_true,Kc_true,...
% lc_true,rr_true,Kr_true,lr_true] = generate_fake_realistic_data(V0, tspan, ...
% sampled_times_individual_data, sampled_times_sum_data, standard_deviation_noise);

%V0 = nan(number_of_different_ratios_tried,number_of_individual_variables);

for iteration = 1:number_of_iterations

    %% 1st step of the ABC algorithm
    % Sampling parameters for the prior
    
    rc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_rc);
    Kc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_Kc);
    lc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_lc);
    rr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_rr);
    Kr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_Kr);
    lr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_lr);
    V0_sum = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_V0_sum);
    
    % Initial values for sum squared errors(iteration_with_certain_ratio)
    
    sum_squared_errors_V_1 = 0;
    sum_squared_errors_V_2 = 0;
    sum_squared_errors_sum = 0;
    
    %% 2nd step of the ABC algorithm
    
    for experiment_with_certain_ratio=1:number_of_different_ratios_tried
        
        sum_of_odds = odds_of_initial_values(experiment_with_certain_ratio,1) + ...
            odds_of_initial_values(experiment_with_certain_ratio,2);
        ratio_V0_C = odds_of_initial_values(experiment_with_certain_ratio,1) / ...
            sum_of_odds;
        ratio_V0_R = odds_of_initial_values(experiment_with_certain_ratio,2) / ...
            sum_of_odds;
        
        V0 = [V0_sum*ratio_V0_C, V0_sum*ratio_V0_R]; 
            

        [t,V_simulated] = ode45(@(t,V) ode_model(t, V, rc, Kc,lc,rr,Kr,lr), ...
        tspan,V0);

        % Visualisation
        % hold on
        % plot(t,V);
        % legend({'Vc','Vr'});
        % hold off


        V_simulated_1 = transpose(V_simulated(:,1)) + standard_deviation_noise*...
            randn(1,length(V_simulated(:,1)));
        V_simulated_2 = transpose(V_simulated(:,2)) + standard_deviation_noise*...
            randn(1,length(V_simulated(:,2)));
        V_simulated_sum = transpose(V_simulated(:,1)) + transpose(V_simulated(:,2)) + ...
            standard_deviation_sum_normal*...
            randn(1,length(V_simulated(:,1) + V_simulated(:,2)));
    
        V_simulated_1 = interp1(t,V_simulated_1,sampled_times_individual_data);
        V_simulated_2 = interp1(t,V_simulated_2,sampled_times_individual_data);
        V_simulated_sum = interp1(t,V_simulated_sum, sampled_times_sum_data);
    
        sum_squared_errors_V_1 = sum_squared_errors_V_1 + ...
            sum((V_simulated_1 - V_individual_control_data_PC3_0Gy(experiment_with_certain_ratio)).^2);
        sum_squared_errors_V_2 = sum_squared_errors_V_2 + ...
            sum((V_simulated_2 - V_individual_resistant_data_PC3_0Gy(experiment_with_certain_ratio)).^2);
        sum_squared_errors_sum = sum_squared_errors_sum + ...
            sum((V_simulated_sum - V_sum_data_average_PC3_0Gy(experiment_with_certain_ratio))...
            .^2);
        
    end
    
    %total_squared_errors = sum_squared_errors_V_1 + sum_squared_errors_V_2 + ...
    %    sum_squared_errors_sum/number_of_individual_variables;
    
    total_squared_errors = sum_squared_errors_V_1 + sum_squared_errors_V_2 + ...
        sum_squared_errors_sum/number_of_individual_variables;
    
    if (total_squared_errors <= tolerance)
        accepted_rc_array(iteration) = rc;
        accepted_Kc_array(iteration) = Kc;
        accepted_lc_array(iteration) = lc;
        accepted_rr_array(iteration) = rr;
        accepted_Kr_array(iteration) = Kr;
        accepted_lr_array(iteration) = lr;
        accepted_V0_sum_array(iteration) = V0_sum;
    end
end

approximate_posterior_mean_rc = nanmean(accepted_rc_array);
approximate_posterior_mean_Kc = nanmean(accepted_Kc_array);
approximate_posterior_mean_lc = nanmean(accepted_lc_array);
approximate_posterior_mean_rr = nanmean(accepted_rr_array);
approximate_posterior_mean_Kr = nanmean(accepted_Kr_array);
approximate_posterior_mean_lr = nanmean(accepted_lr_array);
approximate_posterior_mean_V0_sum = nanmean(accepted_V0_sum_array);

approximate_posterior_median_rc = nanmedian(accepted_rc_array);
approximate_posterior_median_Kc = nanmedian(accepted_Kc_array);
approximate_posterior_median_lc = nanmedian(accepted_lc_array);
approximate_posterior_median_rr = nanmedian(accepted_rr_array);
approximate_posterior_median_Kr = nanmedian(accepted_Kr_array);
approximate_posterior_median_lr = nanmedian(accepted_lr_array);
approximate_posterior_median_V0_sum = nanmedian(accepted_V0_sum_array);

% MAP estimates
MAP_estimate_rc = get_MAP_estimate(accepted_rc_array);
MAP_estimate_Kc = get_MAP_estimate(accepted_Kc_array);
MAP_estimate_lc = get_MAP_estimate(accepted_lc_array);
MAP_estimate_rr = get_MAP_estimate(accepted_rr_array);
MAP_estimate_Kr = get_MAP_estimate(accepted_Kr_array);
MAP_estimate_lr = get_MAP_estimate(accepted_lr_array);
MAP_estimate_V0_sum = get_MAP_estimate(accepted_V0_sum_array);

number_of_accepted_rc = sum(~isnan(accepted_rc_array));
number_of_accepted_Kc = sum(~isnan(accepted_Kc_array));
number_of_accepted_lc = sum(~isnan(accepted_lc_array));
number_of_accepted_rr = sum(~isnan(accepted_rr_array));
number_of_accepted_Kr = sum(~isnan(accepted_Kr_array));
number_of_accepted_lr = sum(~isnan(accepted_lr_array));
number_of_accepted_V0_sum = sum(~isnan(accepted_V0_sum_array));

display('Approximate posterior means for parameters');
display(['rc = ', num2str(approximate_posterior_mean_rc)]);
display(['Kc = ', num2str(approximate_posterior_mean_Kc)]);
display(['lc = ', num2str(approximate_posterior_mean_lc)]);
display(['rr = ', num2str(approximate_posterior_mean_rr)]);
display(['Kr = ', num2str(approximate_posterior_mean_Kr)]);
display(['lr = ', num2str(approximate_posterior_mean_lr)]);
display(['V0_sum = ', num2str(approximate_posterior_mean_V0_sum)]);

display('Approximate posterior medians for parameters');
display(['rc = ', num2str(approximate_posterior_median_rc)]);
display(['Kc = ', num2str(approximate_posterior_median_Kc)]);
display(['lc = ', num2str(approximate_posterior_median_lc)]);
display(['rr = ', num2str(approximate_posterior_median_rr)]);
display(['Kr = ', num2str(approximate_posterior_median_Kr)]);
display(['lr = ', num2str(approximate_posterior_median_lr)]);
display(['V0_sum = ', num2str(approximate_posterior_median_V0_sum)]);

display('Maximum A Posteriori (MAP) estimate for parameters');
display(['rc = ', num2str(MAP_estimate_rc)]);
display(['Kc = ', num2str(MAP_estimate_Kc)]);
display(['lc = ', num2str(MAP_estimate_lc)]);
display(['rr = ', num2str(MAP_estimate_rr)]);
display(['Kr = ', num2str(MAP_estimate_Kr)]);
display(['lr = ', num2str(MAP_estimate_lr)]);
display(['V0_sum = ', num2str(MAP_estimate_V0_sum)]);

display('Number of accepted parameters');
display(['For rc: ', num2str(number_of_accepted_rc)]);
display(['For Kc: ', num2str(number_of_accepted_Kc)]);
display(['For lc: ', num2str(number_of_accepted_lc)]);
display(['For rr: ', num2str(number_of_accepted_rr)]);
display(['For Kr: ', num2str(number_of_accepted_Kr)]);
display(['For lr: ', num2str(number_of_accepted_lr)]);
display(['For V0 sum: ', num2str(number_of_accepted_V0_sum)]);

sum_squared_error_estimated_parameters_mean_posterior = (rc_true_PC3_0Gy - approximate_posterior_mean_rc)^2 + ...
    (Kc_true_PC3_0Gy - approximate_posterior_mean_Kc)^2 + ...
    (lc_true_PC3_0Gy - approximate_posterior_mean_lc)^2 + ...
    (rr_true_PC3_0Gy - approximate_posterior_mean_rr)^2 + ...
    (Kr_true_PC3_0Gy - approximate_posterior_mean_Kr)^2 + ...
    (lr_true_PC3_0Gy - approximate_posterior_mean_lr)^2 + ...
    (V0_sum_true_PC3_0Gy - approximate_posterior_mean_V0_sum)^2;

sum_squared_error_estimated_parameters_median_posterior = (rc_true_PC3_0Gy - approximate_posterior_median_rc)^2 + ...
    (Kc_true_PC3_0Gy - approximate_posterior_median_Kc)^2 + ...
    (lc_true_PC3_0Gy - approximate_posterior_median_lc)^2 + ...
    (rr_true_PC3_0Gy - approximate_posterior_median_rr)^2 + ...
    (Kr_true_PC3_0Gy - approximate_posterior_median_Kr)^2 + ...
    (lr_true_PC3_0Gy - approximate_posterior_median_lr)^2 + ...
    (V0_sum_true_PC3_0Gy - approximate_posterior_median_V0_sum)^2;

sum_squared_error_estimated_parameters_MAP = (rc_true_PC3_0Gy - MAP_estimate_rc)^2 + ...
    (Kc_true_PC3_0Gy - MAP_estimate_Kc)^2 + ...
    (lc_true_PC3_0Gy - MAP_estimate_lc)^2 + ...
    (rr_true_PC3_0Gy - MAP_estimate_rr)^2 + ...
    (Kr_true_PC3_0Gy - MAP_estimate_Kr)^2 + ...
    (lr_true_PC3_0Gy - MAP_estimate_lr)^2 + ...
    (V0_sum_true_PC3_0Gy - MAP_estimate_V0_sum)^2;

display(['sum of squared_error for estimated parameters by mean posterior =', ...
    num2str(sum_squared_error_estimated_parameters_mean_posterior)]);

display(['sum of squared_error for estimated parameters by median posterior =', ...
    num2str(sum_squared_error_estimated_parameters_median_posterior)]);

display(['sum of squared_error for estimated parameters by MAP =', ...
    num2str(sum_squared_error_estimated_parameters_MAP)]);

figure;
subplot(number_of_parameters_to_estimate,1,1);

h = histogram(accepted_rc_array(~isnan(accepted_rc_array)));
morebins(h);
hold on
plot([rc_true_PC3_0Gy,rc_true_PC3_0Gy],[0,10],'r','LineWidth',4);
title('Posterior distribution for rc');
hold off

subplot(number_of_parameters_to_estimate,1,2);
histogram(accepted_Kc_array(~isnan(accepted_Kc_array)));
hold on
plot([Kc_true_PC3_0Gy,Kc_true_PC3_0Gy],[0,10],'r','LineWidth',4);
title('Posterior distribution for Kc');
hold off

subplot(number_of_parameters_to_estimate,1,3);
histogram(accepted_lc_array(~isnan(accepted_lc_array)));
hold on
plot([lc_true_PC3_0Gy,lc_true_PC3_0Gy],[0,10],'r','LineWidth',4);
title('Posterior distribution for lc');
hold off

subplot(number_of_parameters_to_estimate,1,4);
histogram(accepted_rr_array(~isnan(accepted_rr_array)));
hold on
plot([rr_true_PC3_0Gy,rr_true_PC3_0Gy],[0,10],'r','LineWidth',4);
title('Posterior distribution for rr');
hold off

subplot(number_of_parameters_to_estimate,1,5);

histogram(accepted_Kr_array(~isnan(accepted_Kr_array)));
hold on
plot([Kr_true_PC3_0Gy,Kr_true_PC3_0Gy],[0,10],'r','LineWidth',4);
title('Posterior distribution for Kr');
hold off

subplot(number_of_parameters_to_estimate,1,6);

histogram(accepted_lr_array(~isnan(accepted_lr_array)));
hold on
plot([lr_true_PC3_0Gy,lr_true_PC3_0Gy],[0,10],'r','LineWidth',4);
hold off
title('Posterior distribution for lr');

subplot(number_of_parameters_to_estimate,1,7);

histogram(accepted_V0_sum_array(~isnan(accepted_V0_sum_array)));
hold on
plot([V0_sum_true_PC3_0Gy,V0_sum_true_PC3_0Gy],[0,10],'r','LineWidth',4);
hold off
title('Posterior distribution for V0 sum');
end

function parameter = get_positive_parameter_sampled_from_cauchy_distribution(...
    cauchy_prior,sample_size)

    % Default parameter value
    if nargin < 2
        sample_size = 1;
    end
    
    parameter_negative = true;
    
    while parameter_negative
    parameter = random(cauchy_prior,1,sample_size);
    
        if ~(parameter < 0)
            parameter_negative = false;
        end
    end
end

function MAP_estimate = get_MAP_estimate(parameter_array)
[N,edges] = histcounts(parameter_array(~isnan(parameter_array)));
[~,argmax] = max(N);
MAP_estimate = (edges(argmax) + edges(argmax+1))/2;
end
