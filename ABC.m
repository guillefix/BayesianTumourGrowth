function ABC()

number_of_parameters = 6;

%% 1st step of the ABC algorithm

% Boundaries on parameters
lower_bound_rc = 0;
upper_bound_rc = 1;

lower_bound_Kc = 0;
upper_bound_Kc = 1;

lower_bound_lc = 0;
upper_bound_lc = 1;

lower_bound_rr = 0;
upper_bound_rr = 1;

lower_bound_Kr = 0;
upper_bound_Kr = 1;

lower_bound_lr = 0;
upper_bound_lr = 1;

% Sample parameters from uniform prior
rc =  lower_bound_rc + ...
    (upper_bound_rc-lower_bound_rc).*rand;

Kc =  lower_bound_Kc + ...
    (upper_bound_Kc-lower_bound_Kc).*rand;

lc =  lower_bound_lc + ...
    (upper_bound_lc-lower_bound_lc).*rand;

rr =  lower_bound_rr + ...
    (upper_bound_rr-lower_bound_rr).*rand;

Kr =  lower_bound_Kr + ...
    (upper_bound_Kr-lower_bound_Kr).*rand;

lr =  lower_bound_lr + ...
    (upper_bound_lr-lower_bound_lr).*rand;


%% 2nd step of the ABC algorithm

V0 = [1, 1];
tspan = [0 , 100];

[t,V] = ode45(@(t,V) ode_model(t, V, rc, Kc,lc,rr,Kr,lr), ...
tspan,V0);

hold on
plot(t,V);
legend({'Vc','Vr'});
hold off



end
