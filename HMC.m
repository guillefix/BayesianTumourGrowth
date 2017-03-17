x= [0.5,0.5,0.5,0.5,0.5,0.5];
V0=[0.031;0.015];
tspan=[0,15];
sampled_times=1:0.5:15;
standard_deviation_noise=0.01;
[Vdata1, Vdata2] = generate_fake_data(V0, tspan, sampled_times,standard_deviation_noise);
gradE = @(x) gradEdata(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan);
findE = @(x) findEdata(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan);
g = gradE(x);
E = findE(x);
epsilon=0.005;
xs=[];
L=100;
Tau =100;
for l = 1:L
    l
    p = 1*randn(length(x),1); %0.01*
%     pc=p([1,2,4,5]);
%     pc(pc<0) = 0;
%     p([1,2,4,5]) = pc;
    H = p' * p / 2 + E

    xnew = x';
%     xnew([1,2,4,5]) = abs(xnew([1,2,4,5]));
    gnew = g;
    for tau = 1:Tau
        p = p - epsilon * gnew / 2 ;
%         pc=p([1,2,4,5]);
%         pc(pc<0) = 0;
%         p([1,2,4,5]) = pc;
        xnew = xnew + epsilon * p;
%         xnew([1,2,4,5]) = abs(xnew([1,2,4,5]));
        gnew = gradE(xnew);
        p = p - epsilon * gnew / 2 ; 
%         pc=p([1,2,4,5]);
%         pc(pc<0) = 0;
%         p([1,2,4,5]) = pc;
    end
    
    Enew = findE(xnew);
    Hnew = p' * p / 2 + Enew
    dH = Hnew - H ;
    if ( dH <= 0 ) accept = 1;
    elseif (rand() < exp(-dH)) accept = 1;
    else accept = 0;
    end
    if (accept)
        g = gnew;
        x = xnew';
        xs = [xs; x]
        E = Enew;
    end
end
sizes=size(xs);
xs(:,[1,2,4,5]) = xs(:,[1,2,4,5]).*xs(:,[1,2,4,5]);
sum(xs,1)/(sizes(1))