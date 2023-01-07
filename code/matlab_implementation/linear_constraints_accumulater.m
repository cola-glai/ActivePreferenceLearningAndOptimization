
function [A_new,b_new]=linear_constraints_accumulater(x,y,A,b)

    % size(A)
    % size(b)

    [l,d]=size(x);
    h = x.*repmat(y,1,d);
    A_curr = -[ h -h ones(l,1) ];
    b_curr = -ones(l,1);
    % disp('curr:');
    % size(A_curr)
    % size(b_curr)
    
    A_new=[A;A_curr];
    b_new=[b;b_curr];
    
end