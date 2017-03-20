function ABC()

number_of_iterations = 1000;
beginning_time = 0;
end_time = 15;
tspan = [beginning_time , end_time];
sampling_timestep_individual_data = 5;
sampling_timestep_sum_data = 1;
sampled_times_individual_data = beginning_time:sampling_timestep_individual_data:...
    end_time;
sampled_times_sum_data = beginning_time:sampling_timestep_sum_data:end_time;
number_of_individual_variables = 2;
number_of_sampled_times_individual_variables = size(sampled_times_individual_data,2);
number_of_sampled_times_sum = size(sampled_times_sum_data,2);
total_number_of_sampled_times = number_of_individual_variables *...
    number_of_sampled_times_individual_variables + ...
    number_of_sampled_times_sum;
tolerance = 0.1*total_number_of_sampled_times;
standard_deviation_noise = 0.01;

accepted_rc_array = nan(number_of_iterations,1);
accepted_Kc_array = nan(number_of_iterations,1);
accepted_lc_array = nan(number_of_iterations,1);
accepted_rr_array = nan(number_of_iterations,1);
accepted_Kr_array = nan(number_of_iterations,1);
accepted_lr_array = nan(number_of_iterations,1);

V0 = [1, 2];

% Boundaries on parameters
lower_bound_rc = 0;
upper_bound_rc = 3;

lower_bound_Kc = 0;
upper_bound_Kc = 3;

lower_bound_lc = 0;
upper_bound_lc = 3;

lower_bound_rr = 0;
upper_bound_rr = 3;

lower_bound_Kr = 0;
upper_bound_Kr = 3;

lower_bound_lr = 0;
upper_bound_lr = 3;



for iteration = 1:number_of_iterations

    %% 1st step of the ABC algorithm

    % Sample parameters from uniform prior
    rc =  lower_bound_rc + ...
        (upper_bound_rc-lower_bound_rc).*rand();

    Kc =  lower_bound_Kc + ...
        (upper_bound_Kc-lower_bound_Kc).*rand();
    
    lc =  lower_bound_lc + ...
        (upper_bound_lc-lower_bound_lc).*rand();

    rr =  lower_bound_rr + ...
        (upper_bound_rr-lower_bound_rr).*rand();

    Kr =  lower_bound_Kr + ...
        (upper_bound_Kr-lower_bound_Kr).*rand();

    lr =  lower_bound_lr + ...
        (upper_bound_lr-lower_bound_lr).*rand();


    %% 2nd step of the ABC algorithm

    [t,V_simulated] = ode45(@(t,V) ode_model(t, V, rc, Kc,lc,rr,Kr,lr), ...
    tspan,V0);

    % Visualisationotal_squared_errors
    % hold on
    % plot(t,V);
    % legend({'Vc','Vr'});
    % hold off

    [V_fakedata_sum,V_fakedata_1, V_fakedata_2] = generate_fake_realistic_data(V0, tspan, ...
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

display('Number of accepted parameters');
display(['For rc: ', num2str(number_of_accepted_rc)]);
display(['For Kc: ', num2str(number_of_accepted_Kc)]);
display(['For lc: ', num2str(number_of_accepted_lc)]);
display(['For rr: ', num2str(number_of_accepted_rr)]);
display(['For Kr: ', num2str(number_of_accepted_Kr)]);
display(['For lr: ', num2str(number_of_accepted_lr)]);

end

