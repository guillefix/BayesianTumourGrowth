rc=1;
Kc=1;
lc=1;
rr=2;
Kr=2;
lr=2;

dvdt = @(t,V) [V(1)*(rc*(1-V(1)/Kc)-lr*V(2));V(2)*(rr*(1-V(2)/Kr)-lc*V(1))] ;
V0=[1;1];
tspan=[0,5];
[t,V]=ode45(dvdt,tspan,V0);
plot(t,V(:,1))
hold on
plot(t,V(:,2))

Vdata1=V(1:5:end,1)
Vdata2=V(1:5:end,2)
t(1:5:end)