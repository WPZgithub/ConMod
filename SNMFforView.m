function [ H ] = SNMFforView(X, Hc, H, lambda, maxIter, epsilon)
% Multi-View Non-negative symmetric Matrix Factorization for each view
% 
% INPUT:
%   X: the adjacency matrix of a network
%   Hc: initialization for consensus factor matrix 
%   H: initialization for factor matrix of each view
%   lambda: a vector containing the parameters for balancing the relative
%           weight among different views
%   MaxIter: the maximal number of iterations for alternating minimization
%   epsilon: the convergence parameter
%
% OUTPUT:
%   H: the factor matrix
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

obj_old = norm(X - H*H', 'fro')^2 + lambda*norm(H - Hc, 'fro')^2;

% Fixing Hc, minimize objective function over H
for iter = 1:maxIter   
    % Update rule
    temp_1 = 2*X*H + lambda*Hc;
    temp_2 = 2*H*(H'*H) + lambda*H;
    H = H.*(temp_1./(temp_2+eps));
    
    % Objective function
    obj_body = norm(X - H*H', 'fro')^2;
    obj_consensus = norm(H - Hc, 'fro')^2;
    obj = obj_body + lambda*obj_consensus;

    Delta = obj_old - obj;
    obj_old = obj;
    
    if Delta < epsilon
        break;
    end
end

end