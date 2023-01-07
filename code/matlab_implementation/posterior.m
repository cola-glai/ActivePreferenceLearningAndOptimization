clear all;
clc;
close all;

X = load('pop_zdt1.txt');
% n -> number of objects
% d -> dimension
[n, d] = size(X);
w_star = [0.382039 0.382046];

active_ranking(X, w_star);

% fprintf('<info> w_star: ');
% disp(w_star);
% fprintf('<info> X:\n');
% disp(X);
