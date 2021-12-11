%% MU-ANALYSIS
clear all; close all;

%System
A = [0 1; 2 -1];

yalmip('clear');
Const = [];
eta = .0001;
n = size(A,1);
X = sdpvar(n,n);
gamma2 = sdpvar(1);
Const = [Const; X>=eta*eye(size(X))];
mat1 = A'*X*A-gamma2*X;
Const = [Const; mat1 <= -eta*eye(size(mat1))];

opt=sdpsettings('solver','sedumi','verbose',0);
optimize(Const, gamma2, opt);
XX = value(X);
gamma = sqrt(value(gamma2));

disp(XX);
disp(gamma);