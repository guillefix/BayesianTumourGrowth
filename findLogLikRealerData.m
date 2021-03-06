function E = findLogLikRealerData(x,odds,V_data_sums,sigmaDataSums,Vdata1s,sigmaData1s,Vdata2s,sigmaData2s,N,V0_sum,sampled_times_individual_data, ...
    sampled_times_sum_data, standard_deviation_noise,tspan)

    E = 0;

    for ii=1:length(odds)
        o=odds(ii,:);
        V0=o*V0_sum/sum(o);

        sigmaDataSum = sigmaDataSums(ii,:);
        sigmaData1 = sigmaData1s(ii,:);
        sigmaData2 = sigmaData2s(ii,:);

        [t,V]=ode45(@(t,V) ode_model(t,V, x(1),x(2),x(3),x(4),x(5),x(6)),tspan,V0);

        V1=interp1(t,V(:,1),sampled_times_individual_data);
        V2=interp1(t,V(:,2),sampled_times_individual_data);

        V_sum = interp1(t,V(:,1) + V(:,2),sampled_times_sum_data);

        sigma = standard_deviation_noise/sqrt(N);

        T=length(sampled_times_individual_data);

        for i=1:T
            if (i==T) sigma = standard_deviation_noise/sqrt(N/2);
            else sigma = standard_deviation_noise/sqrt(N);
            end
            E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-((Vdata1s(ii,i)-V1(i))^2+sigmaData1(i)^2)/(2*sigma^2)));
            E = E - log(1/(sqrt(2*pi*sigma^2))*exp(-((Vdata2s(ii,i)-V2(i))^2+sigmaData2(i)^2)/(2*sigma^2)));
        end

        for i=1:length(sampled_times_sum_data)
            if (i==T) sigma = standard_deviation_noise/sqrt(N/2);
            else sigma = standard_deviation_noise/sqrt(N);
            end
            E = E - log(1/(sqrt(2*pi*2*sigma^2))*exp(-((V_data_sums(ii,i)-V_sum(i))^2+sigmaDataSum(i)^2)/(2*2*sigma^2))); % variance is twice for sum
        end

    end

end
