clear all;
clc;
close all;

X = load('pop.txt');
% n -> number of objects
% d -> dimension
[n, d] = size(X);
w_star = [0.382 0.382];

[w_approximate]=active_ranking(X, w_star);
fid=fopen('w_approximate.txt', 'w');
fprintf(fid, '%f %f\n', w_approximate(1), w_approximate(2));
fclose(fid);

% fprintf('<info> w_star: ');
% disp(w_star);
% fprintf('<info> X:\n');
% disp(X);
