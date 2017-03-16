x= [1,1,1,2,2,2,1];
[Vdata1, Vdata2] = generate_fake_data();
gradE = @(x) gradEdata(x,Vdata1,Vdata2);
findE = @(x) findEdata(x,Vdata1,Vdata2);
g = gradE(x);
E = findE(x);
epsilon=0.055;
xs=[];
L=50;
Tau =19;
for l = 1:L
    l
    p = 0.01*randn(length(x),1);
    H = p' * p / 2 + E

    xnew = x';
    gnew = g;
    for tau = 1:Tau
        p = p - epsilon * gnew / 2 ; 
        xnew = xnew + epsilon * p;
        gnew = gradE(xnew);
        p = p - epsilon * gnew / 2 ; 
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
sum(xs,1)/(sizes(1))