function g = gradEdata(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan)

h=0.0005;
g=zeros(length(x),1);
for i=1:length(x)
    xh = x;
    xh(i) = xh(i) +h;
    g(i) = (findEdata(xh,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan)...
        -findEdata(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan))/h;
end

end