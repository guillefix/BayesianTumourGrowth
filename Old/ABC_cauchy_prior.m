function ABC()

number_of_iterations = 3000;
beginning_time = 0;
end_time = 15;
number_of_parameters = 6;
number_of_individual_variables = 2;

tspan = [beginning_time , end_time];
sampling_timestep_individual_data = 5;
sampling_timestep_sum_data = 1;
sampled_times_individual_data = beginning_time:sampling_timestep_individual_data:...
    end_time;
sampled_times_sum_data = beginning_time:sampling_timestep_sum_data:end_time;


number_of_sampled_times_individual_variables = size(sampled_times_individual_data,2);
number_of_sampled_times_sum = size(sampled_times_sum_data,2);

total_number_of_sampled_times = number_of_individual_variables *...
    number_of_sampled_times_individual_variables + ...
    number_of_sampled_times_sum;

tolerance = 0.01*total_number_of_sampled_times;
standard_deviation_noise = 0.01;

accepted_rc_array = nan(number_of_iterations,1);
accepted_Kc_array = nan(number_of_iterations,1);
accepted_lc_array = nan(number_of_iterations,1);
accepted_rr_array = nan(number_of_iterations,1);
accepted_Kr_array = nan(number_of_iterations,1);
accepted_lr_array = nan(number_of_iterations,1);

V0 = [1, 2];

% Generating Cauchy prior
cauchy_prior_rc = makedist('tLocationScale','mu',0.236,'sigma',1,'nu',1);
cauchy_prior_Kc = makedist('tLocationScale','mu',0.473,'sigma',1,'nu',1);
cauchy_prior_lc = makedist('tLocationScale','mu',0.655,'sigma',1,'nu',1);
cauchy_prior_rr = makedist('tLocationScale','mu',0.401,'sigma',1,'nu',1);
cauchy_prior_Kr = makedist('tLocationScale','mu',0.773,'sigma',1,'nu',1);
cauchy_prior_lr = makedist('tLocationScale','mu',0.480,'sigma',1,'nu',1);

% % Parameters for normal prior
% mean_rc = 2;
% standard_deviation_rc = 1;
% 
% mean_Kc = 2;
% standard_deviation_Kc = 1;
% 
% mean_lc = 2;
% standard_deviation_lc = 1;
% 
% mean_rr = 2;
% standard_deviation_rr = 1;
% 
% mean_Kr = 2;
% standard_deviation_Kr = 1;
% 
% mean_lr = 2;
% standard_deviation_lr = 1;



for iteration = 1:number_of_iterations

    %% 1st step of the ABC algorithm
    % Sampling parameters for the prior
    
    rc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_rc);
    Kc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_Kc);
    lc = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_lc);
    rr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_rr);
    Kr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_Kr);
    lr = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior_lr);
    
    %% 2nd step of the ABC algorithm

    [t,V_simulated] = ode45(@(t,V) ode_model(t, V, rc, Kc,lc,rr,Kr,lr), ...
    tspan,V0);

    % Visualisationotal_squared_errors
    % hold on
    % plot(t,V);
    % legend({'Vc','Vr'});
    % hold off

    [V_fakedata_sum,V_fakedata_1, V_fakedata_2,rc_true,Kc_true,...
    lc_true,rr_true,Kr_true,lr_true] = generate_fake_realistic_data(V0, tspan, ...
    sampled_times_individual_data, sampled_times_sum_data, standard_deviation_noise);
    
    V_simulated_1_no_noise = interp1(t,V_simulated(:,1),sampled_times_individual_data);
    V_simulated_2_no_noise = interp1(t,V_simulated(:,2),sampled_times_individual_data);
    V_simulated_sum_no_noise = interp1(t,V_simulated(:,1) + V_simulated(:,2), ...
        sampled_times_sum_data);
    
    V_simulated_1 = V_simulated_1_no_noise + standard_deviation_noise*...
        randn(1,number_of_sampled_times_individual_variables);
    V_simulated_2 = V_simulated_2_no_noise + standard_deviation_noise*...
        randn(1,number_of_sampled_times_individual_variables);
    V_simulated_sum = V_simulated_sum_no_noise + standard_deviation_noise*...
        randn(1,number_of_sampled_times_sum);
    
    sum_squared_errors_V_1 = sum((V_simulated_1 - V_fakedata_1).^2);
    sum_squared_errors_V_2 = sum((V_simulated_2 - V_fakedata_2).^2);
    sum_squared_errors_sum = sum((V_simulated_sum - V_fakedata_sum).^2);
    
    total_squared_errors = sum_squared_errors_V_1 + sum_squared_errors_V_2 + ...
        sum_squared_errors_sum;
    
    if (total_squared_errors <= tolerance)
        accepted_rc_array(iteration) = rc;
        accepted_Kc_array(iteration) = Kc;
        accepted_lc_array(iteration) = lc;
        accepted_rr_array(iteration) = rr;
        accepted_Kr_array(iteration) = Kr;
        accepted_lr_array(iteration) = lr;
    end
end

approximate_posterior_mean_rc = nanmean(accepted_rc_array);
approximate_posterior_mean_Kc = nanmean(accepted_Kc_array);
approximate_posterior_mean_lc = nanmean(accepted_lc_array);
approximate_posterior_mean_rr = nanmean(accepted_rr_array);
approximate_posterior_mean_Kr = nanmean(accepted_Kr_array);
approximate_posterior_mean_lr = nanmean(accepted_lr_array);

approximate_posterior_median_rc = nanmedian(accepted_rc_array);
approximate_posterior_median_Kc = nanmedian(accepted_Kc_array);
approximate_posterior_median_lc = nanmedian(accepted_lc_array);
approximate_posterior_median_rr = nanmedian(accepted_rr_array);
approximate_posterior_median_Kr = nanmedian(accepted_Kr_array);
approximate_posterior_median_lr = nanmedian(accepted_lr_array);

% MAP estimates
MAP_estimate_rc = get_MAP_estimate(accepted_rc_array);
MAP_estimate_Kc = get_MAP_estimate(accepted_Kc_array);
MAP_estimate_lc = get_MAP_estimate(accepted_lc_array);
MAP_estimate_rr = get_MAP_estimate(accepted_rr_array);
MAP_estimate_Kr = get_MAP_estimate(accepted_Kr_array);
MAP_estimate_lr = get_MAP_estimate(accepted_lr_array);

number_of_accepted_rc = sum(~isnan(accepted_rc_array));
number_of_accepted_Kc = sum(~isnan(accepted_Kc_array));
number_of_accepted_lc = sum(~isnan(accepted_lc_array));
number_of_accepted_rr = sum(~isnan(accepted_rr_array));
number_of_accepted_Kr = sum(~isnan(accepted_Kr_array));
number_of_accepted_lr = sum(~isnan(accepted_lr_array));

display('Approximate posterior means for parameters');
display(['rc = ', num2str(approximate_posterior_mean_rc)]);
display(['Kc = ', num2str(approximate_posterior_mean_Kc)]);
display(['lc = ', num2str(approximate_posterior_mean_lc)]);
display(['rr = ', num2str(approximate_posterior_mean_rr)]);
display(['Kr = ', num2str(approximate_posterior_mean_Kr)]);
display(['lr = ', num2str(approximate_posterior_mean_lr)]);

display('Approximate posterior medians for parameters');
display(['rc = ', num2str(approximate_posterior_median_rc)]);
display(['Kc = ', num2str(approximate_posterior_median_Kc)]);
display(['lc = ', num2str(approximate_posterior_median_lc)]);
display(['rr = ', num2str(approximate_posterior_median_rr)]);
display(['Kr = ', num2str(approximate_posterior_median_Kr)]);
display(['lr = ', num2str(approximate_posterior_median_lr)]);

display('Maximum A Posteriori (MAP) estimate for parameters');
display(['rc = ', num2str(MAP_estimate_rc)]);
display(['Kc = ', num2str(MAP_estimate_Kc)]);
display(['lc = ', num2str(MAP_estimate_lc)]);
display(['rr = ', num2str(MAP_estimate_rr)]);
display(['Kr = ', num2str(MAP_estimate_Kr)]);
display(['lr = ', num2str(MAP_estimate_lr)]);

display('Number of accepted parameters');
display(['For rc: ', num2str(number_of_accepted_rc)]);
display(['For Kc: ', num2str(number_of_accepted_Kc)]);
display(['For lc: ', num2str(number_of_accepted_lc)]);
display(['For rr: ', num2str(number_of_accepted_rr)]);
display(['For Kr: ', num2str(number_of_accepted_Kr)]);
display(['For lr: ', num2str(number_of_accepted_lr)]);

sum_squared_error_estimated_parameters_mean_posterior = (rc_true - approximate_posterior_mean_rc)^2 + ...
    (Kc_true - approximate_posterior_mean_Kc)^2 + ...
    (lc_true - approximate_posterior_mean_lc)^2 + ...
    (rr_true - approximate_posterior_mean_rr)^2 + ...
    (Kr_true - approximate_posterior_mean_Kr)^2 + ...
    (lr_true - approximate_posterior_mean_lr)^2;

sum_squared_error_estimated_parameters_median_posterior = (rc_true - approximate_posterior_median_rc)^2 + ...
    (Kc_true - approximate_posterior_median_Kc)^2 + ...
    (lc_true - approximate_posterior_median_lc)^2 + ...
    (rr_true - approximate_posterior_median_rr)^2 + ...
    (Kr_true - approximate_posterior_median_Kr)^2 + ...
    (lr_true - approximate_posterior_median_lr)^2;

sum_squared_error_estimated_parameters_MAP = (rc_true - MAP_estimate_rc)^2 + ...
    (Kc_true - MAP_estimate_Kc)^2 + ...
    (lc_true - MAP_estimate_lc)^2 + ...
    (rr_true - MAP_estimate_rr)^2 + ...
    (Kr_true - MAP_estimate_Kr)^2 + ...
    (lr_true - MAP_estimate_lr)^2;
parameter = random(cauchy_prior,1,1);
display(['sum of squared_error for estimated parameters by mean posterior =', ...
    num2str(sum_squared_error_estimated_parameters_mean_posterior)]);

display(['sum of squared_error for estimated parameters by median posterior =', ...
    num2str(sum_squared_error_estimated_parameters_median_posterior)]);

display(['sum of squared_error for estimated parameters by MAP =', ...
    num2str(sum_squared_error_estimated_parameters_MAP)]);

figure;
subplot(6,1,1);

histogram(accepted_rc_array(~isnan(accepted_rc_array)));
hold on
plot([rc_true,rc_true],[0,10],'r','LineWidth',4);
title('Posterior distribution for rc');
hold off

subplot(6,1,2);
histogram(accepted_Kc_array(~isnan(accepted_Kc_array)));
hold on
plot([Kc_true,Kc_true],[0,10],'r','LineWidth',4);
title('Posterior distribution for Kc');
hold off

subplot(6,1,3);
histogram(accepted_lc_array(~isnan(accepted_lc_array)));
hold on
plot([lc_true,lc_true],[0,10],'r','LineWidth',4);
title('Posterior distribution for lc');
hold off

subplot(6,1,4);
histogram(accepted_rr_array(~isnan(accepted_rr_array)));
hold on
plot([rr_true,rr_true],[0,10],'r','LineWidth',4);
title('Posterior distribution for rr');
hold off

subplot(6,1,5);

histogram(accepted_Kr_array(~isnan(accepted_Kr_array)));
hold on
plot([Kr_true,Kr_true],[0,10],'r','LineWidth',4);
title('Posterior distribution for Kr');
hold off

subplot(6,1,6);

histogram(accepted_lr_array(~isnan(accepted_lr_array)));
hold on
plot([lr_true,lr_true],[0,10],'r','LineWidth',4);
hold off
title('Posterior distribution for lr');
end

function parameter = get_positive_parameter_sampled_from_cauchy_distribution(cauchy_prior)
    parameter_negative = true;
    
    while parameter_negative
    parameter = random(cauchy_prior,1,1);
    
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
