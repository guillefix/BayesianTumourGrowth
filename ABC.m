function ABC()

% rc,Kc,lc,rr,Kr,lr

number_of_parameters = 6;
lower_bound_parameters = zeros(1,number_of_parameters);
upper_bound_parameters = ones(1, number_of_parameters);

rc =  lower_bound_parameters + ...
    (upper_bound_parameters-lower_bound_parameters).*rand;

[rc,Kc,lc,rr,Kr,lr] = fun();

V0 = [1, 1];
tspan = [0 , 100];

[t,V] = ode45(@(t,V) ode_model(t, V, rc, Kc,lc,rr,Kr,lr), ...
tspan,V0);




end

