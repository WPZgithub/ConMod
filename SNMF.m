function [ H,objj ] = SNMF(X, K, H_init, maxIter, epsilon)
% Symmetric Non-Negtive Matrix Factorization

% INPUT:
%       X: the adjacency matrix of a network
%       K: the number of hidden factors
%       maxIter: the maximal number of iterations for alternating minimization
%       epsilon: the convergence parameter
%
% OUTPUT:
%       H: the factor matirx
%       objj: the value of objective function
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

%% Normalize the network
X = X/sqrt(trace(X'*X));

%% Initializaiton
N = size(X,1);
if isempty(H_init)
    H = rand(N,K);
else
    H = H_init;
end
H = H/sqrt(trace(H'*H));

obj_old = norm(X - H*H', 'fro')^2;
beta = 0.5;
objj = obj_old;
%% Alternating update
for iter = 1:maxIter  
    temp_1 = X*H;
    temp_2 = (H')*H;
    temp_2 = H*temp_2;
    H = H.*(1 - beta + beta*(temp_1./(temp_2+eps)));
    
    obj = norm(X - H*H', 'fro')^2;
    Delta = obj_old - obj;
    obj_old = obj;
    
    objj = [objj, obj];
    if Delta < epsilon
        break;
    end
end

end