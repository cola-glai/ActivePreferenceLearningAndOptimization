function z=ZDT1(x)

    n=numel(x);
    
    f1=x(1);
    i=2:n;
    g=1+9*sum(x(i))/(n-1);
    f2=1-sqrt(f1/g);
    
    z=[f1 f2]';

end