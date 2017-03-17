function E = findEdata(x,Vdata1, Vdata2,V0,sampled_times,standard_deviation_noise,tspan) 

    [t,V]=ode45(@(t,V) ode_model(t,V, x(1)^2,x(2)^2,x(3),x(4)^2,x(5)^2,x(6)),tspan,V0);

    V1=interp1(t,V(:,1),sampled_times);
    V2=interp1(t,V(:,2),sampled_times);

    sigma = standard_deviation_noise; %1
    E = 0;
    for i=1:length(sampled_times)
%         E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-(Vdata1(i)-V1(i))^2/(2*sigma^2)));
%         E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-(Vdata2(i)-V2(i))^2/(2*sigma^2)));
        E = E+(Vdata1(i)-V1(i))^2/(2*sigma^2);
        E = E+(Vdata2(i)-V2(i))^2/(2*sigma^2);
    end
    
%     if (x(2) < 0 || x(5) < 0 || x(1) < 0 || x(4) < 0)
%         xneg=x([2,5,1,4]);
%         xneg = neg(x<0);
%         E=E+100*norm(xneg)^2;
%     end
    
%     1/(sqrt(2*pi*sigma^2))

end