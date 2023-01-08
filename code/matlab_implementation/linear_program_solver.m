function [w,fval,flag,acc]=linear_program_solver(x,y,append)
% x: (n-x-d) in \R^p
% y: (n-x-1) in {-1,1}
% append: 1-x-1 \in {0,1} indicating if a 1 should be appended to x

if nargin>2 && append
    x = [x ones(size(x,1),1)];
end

x
y

% add by gylai
% d == [origin dimension] + 1
[l,d]=size(x);

% add by gylai
% repmat(y,1,d) -> d cols of y
h = x.*repmat(y,1,d);

A = -[ h -h ones(l,1) ];
b = -ones(l,1);
f = [zeros(2*d,1);1];
LB = [ zeros(2*d,1); -1];
UB = inf(2*d+1,1);

% A = -[ h(:,1:d) -h(:,1:d) ones(l,1) ];
% b = h(:,end);
% f = [zeros(2*d,1);1];
% LB = [ zeros(2*d,1); -1];
% UB = inf(2*d+1,1);

% A = - [ h(:,1) h(:,2:end-1) -h(:,2:end-1) ones(l,1) ];
% b = h(:,end);
% f = [zeros(2*d+1,1);1];
% LB = [ zeros(2*d+1,1); -10];
% UB = inf(2*d+2,1);

if nargin > 2
% problem.x0 = [ w0(1); w0(2:end); -w0(2:end); 1];
end

problem.f = f;
problem.Aineq = A;
problem.bineq = b;
problem.lb = LB;
problem.ub = UB;
problem.solver = 'linprog';
problem.options = optimset('display','off','tolfun',1e-8);
% add by gylai
% u -> optimal solution, the size of u is the same as the size of f
% fval -> fval = f ' * u
% flag -> exit flag
[u, fval, flag ] = linprog(problem);

if isempty(u)
    disp('set tolfun as 1e-1');
    problem.options = optimset('display','off','tolfun',1e-1);
    [u, fval, flag ] = linprog(problem);
end
% fprintf('======= size(u) =========\n');
% disp(size(u));
% add by gylai
% FIXME
u
w = u(1:d)-u(d+1:2*d);
% w = [ u(1); u(2:d+1)-u(d+2:2*d+1)]; 
% add by gylai
% sign(num)
% if num < 0, then sign(num) == -1
% if num > 0, then sign(num) == 1
% if num = 0, then sign(num) == 0
% x -> the hyperplane with (d + 1) parameters
% x*w;
% y;
acc = sum( sign(x*w)==y )/length(y)*100;
w
acc

% uncommment by gylai
% disp(['Solution=' num2str(fval,'%1.3g') '     Accuracy=' num2str(acc,'%1.3g') '     Training data Size=' num2str(l) '     Time=' num2str(time)])
% 


% X = LINPROG(f,A,b) attempts to solve the linear programming problem:
%          
%              min f'*x    subject to:   A*x <= b 
%               x
%  
%     X = LINPROG(f,A,b,Aeq,beq) solves the problem above while additionally
%     satisfying the equality constraints Aeq*x = beq.

%       1  LINPROG converged to a solution X.
%       0  Maximum number of iterations reached.
%      -2  No feasible point found.
%      -3  Problem is unbounded.
%      -4  NaN value encountered during execution of algorithm.
%      -5  Both primal and dual problems are infeasible.
%      -7  Magnitude of search direction became too small; no further 
%           progress can be made. The problem is ill-posed or badly 
%           conditioned.