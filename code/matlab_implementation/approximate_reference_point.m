
function w_approximate=approximate_reference_point(A,b,f,LB,UB,d)
    [nPair,~]=size(A);
    % disp('number of accumulated active pairwise comparisons by the DM:');
    % disp(nPair);
    problem.f = f;
    problem.Aineq = A;
    problem.bineq = b;
    problem.lb = LB;
    problem.ub = UB;
    problem.solver = 'linprog';
    problem.options = optimset('display','off','tolfun',1e-8);
    [u,~,~] = linprog(problem);
    
    if isempty(u)
    problem.options = optimset('display','off','tolfun',1e-1);
    [u,~,~] = linprog(problem);
    end
    
    w = u(1:d)-u(d+1:2*d);
    w_approximate = w(1 : end - 1) / w(end);
    w_approximate = w_approximate';
    
end