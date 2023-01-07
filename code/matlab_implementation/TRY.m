
%% test function compare_sort
% goodInds = [1   1   1   0   0   0   1   0];
% Qhyp = [[0     1     1    -1];
%             [-1     0    -1    -1];
%             [-1     1     0    -1];
%             [1     1     1     0]];
%  index = compare_sort(1:sum(goodInds),Qhyp);

%% test function linear_program_solver
X = [0.5617    0.5006   -0.6587;
    0.5558    0.5616   -0.6129;
   -0.5634   -0.3678    0.7398;
   -0.8803    0.2895    0.3758;
   -0.2467    0.9161   -0.3160;
   -0.9206    0.3607    0.1498]
Y = [1;
     1;
    -1;
    -1;
    -1;
     1]
[w,fval,flag,accuracy_n]=linear_program_solver(X,Y,0);