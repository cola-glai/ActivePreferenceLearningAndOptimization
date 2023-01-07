syms x;
% golden point for ZDT1 when f1 == f2
res=solve(x+sqrt(x)-1,x);
res=vpa(res,6);

% golden point for DTLZ1 when f1 == f2 ==f3
res=solve(3*x-0.5,x);
res=vpa(res,6)