goodInds = [1   1   1   0   0   0   1   0];
Qhyp = [[0     1     1    -1];
            [-1     0    -1    -1];
            [-1     1     0    -1];
            [1     1     1     0]];
 index = compare_sort(1:sum(goodInds),Qhyp)
