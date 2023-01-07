%
% created by gylai on Jan 5, 2022
%
clc;
clear;
close all;
best_solutions=cell(1,5);

tic
parfor i=1:5
    best_solutions{i}=parallel_run();
    disp(best_solutions{i});
end
toc

