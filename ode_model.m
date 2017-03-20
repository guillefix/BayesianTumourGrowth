function dVdt = ode_model(t, V, rc, Kc,lc,rr,Kr,lr)

dVdt = zeros(2,1);

dVdt(1) = V(1)*(rc*(1-V(1)/Kc)-lr*V(2)); % V_c
dVdt(2) = V(2)*(rr*(1-V(2)/Kr)-lc*V(1)); % V_r


