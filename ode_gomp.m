function dVdt = ode_gomp(t, V, rc, Kc,lc,rr,Kr,lr)

dVdt = zeros(2,1);

dVdt(1) = V(1)*(rc*log(Kc/V(1))-lr*V(2)); % V_c
dVdt(2) = V(2)*(rr*log(Kr/V(2))-lc*V(1)); % V_r
end
