function E = findEdata(x,Vdata1, Vdata2) 

if (x(2) < 0 || x(5) < 0 || x(1) < 0 || x(4) < 0) E = 99999999;
else
    V0=[1;2];
    tspan=[0,5];
    [t,V]=ode45(@(t,V) ode_model(t,V, x(1),x(2),x(3),x(4),x(5),x(6)),tspan,V0);

    ts=1:0.01:5;
    V1=interp1(t,V(:,1),ts);
    V2=interp1(t,V(:,2),ts);

    sigma = x(7);
    E = 0;
    for i=1:length(ts)
        E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-(Vdata1(i)-V1(i))^2/(2*sigma^2)));
        E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-(Vdata2(i)-V2(i))^2/(2*sigma^2)));
    end
    
%     1/(sqrt(2*pi*sigma^2))
end

end