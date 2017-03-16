function [dvdt,y] = ode(t, V, rc,Kc,lc,rr,Kr,lr varargin)

dvdt = @(t,V) [V(1)*(rc*(1-V(1)/Kc)-lr*V(2));V(2)*(rr*(1-V(2)/Kr)-lc*V(1))] ;
y=V;

