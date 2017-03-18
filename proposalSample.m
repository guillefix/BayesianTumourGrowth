function xnew = proposalSample(x,sigma) 

    xnew = x + sigma*randn(size(x));

end