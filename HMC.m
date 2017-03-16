g = gradE(x) ;
E = findE(x) ;
findE = @(x) 
epsilon=0.1;
for l = 1:L
    p = randn(size(x)) ;
    is Normal(0,1)
    H = p' * p / 2 + E ;

    xnew = x ;  gnew = g ;
    for tau = 1:Tau
        p = p - epsilon * gnew / 2 ; 
        xnew = xnew + epsilon * p ;  
        gnew = gradE(xnew) ;
        p = p - epsilon * gnew / 2 ; 
    end
    Enew = findE(xnew);
    Hnew = p' * p / 2 + Enew;
    dH = Hnew - H ;
    if ( dH < 0 ) accept = 1;
    elseif (rand() < exp(-dH)) accept = 1;
    else accept = 0;
    end
    if (accept)
        g = gnew;
        x = xnew;
        E = Enew;
    end
end