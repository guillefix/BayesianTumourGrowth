function g = gradEdata(x,Vdata1,Vdata2)

h=0.005;
g=zeros(length(x),1);
for i=1:length(x)
    xh = x;
    xh(i) = xh(i) +h;
    g(i) = (findEdata(xh,Vdata1,Vdata2)-findEdata(x,Vdata1,Vdata2))/h;
end

end