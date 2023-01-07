clear all;
clc;
close all;

% number of objects
n = 100;
% dimension
d = 2;

X = rand(n, d) * 100;
X = [0.1598   0.4231;
   0.9921   0.6555;
   0.0396   0.9304;
   0.5975   0.3444;
   0.5423   0.2546;
   0.0572   0.4447;
   0.6315   0.7902;
   0.4236   0.6954;
   0.8423   0.5577;
   0.9063   0.2874
];
w_star = rand(1, d) *100;
w_star = [0.1458   0.5793];

active_ranking(X, w_star);

fprintf('<info> w_star: ');
disp(w_star);
fprintf('<info> X:\n');
disp(X);
