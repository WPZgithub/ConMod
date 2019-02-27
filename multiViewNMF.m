function [ H, Hc, objValue ] = multiViewNMF( X, K, lambda, maxIter )
% Multi-View Non-negative symmetric Matrix Factorization
%
% INPUT:
%   X: a cell which contains symmetric matrices
%   K: the number of hidden factors
%   lambda: a vector which contains the parameters for balancing the relative
%                weight among different views
%   maxiter: the maximum number of iterations
%
% OUTPUT:
%   H: a cell which contains factor matrices for all views
%   Hc: the result consensus factor matrix
%   objValue: the value of objective function
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

%% Test the legitimacy of the input
Xcount = length(X);
Ncount = size(X{1}, 1);
n = Ncount;
for i = 1:Xcount
    n_old = n;
    [n,m] = size(X{i});
    if min(min(X{i})) < 0
        error 'Input matrix elements can not be negative';
    end
    if n ~= m
        error 'Input matrices should be symmetric';
    end
    if n ~= n_old
        error 'Input matrices should have the same rows';
    end
end

H = cell(Xcount, 1);
Hc = zeros(Ncount, K);

%% Normalize input matrix X{i}
for i = 1:Xcount
    X{i} = X{i}/sqrt(trace(X{i}'*X{i}));
end

%% Initialize H and Hc
iterNum = 30;
if Ncount > 2000
    iterNum = 20;
end
H_init = [];
H{1} = SNMF(X{1}, K, H_init, iterNum, 1e-5);
for i = 2:Xcount
    H{i} = SNMF(X{i}, K, H{i-1}, iterNum, 1e-5);
end
for i = 1:Xcount
    Hc = Hc + lambda(i)*H{i};
end
Hc = Hc/sum(lambda);

obj_old = 0;
for i = 1:Xcount
    obj_body = norm(X{i} - H{i}*H{i}', 'fro')^2;
    obj_consensus = norm(H{i} - Hc, 'fro')^2;
    obj_old = obj_old + obj_body + lambda(i)*obj_consensus;
end

%%  Update process
objValue = obj_old;
for iter = 1:maxIter
    % Fixing Hc, minimize objective function over H
    maxIterforView = 40;
    for i = 1:Xcount
        [H{i}] = SNMFforView(X{i}, Hc, H{i}, lambda(i), maxIterforView, 1e-6);
    end 
    
    % Fixing H, minimize objective function over Hc
    Hc = zeros(Ncount, K);
    for i = 1:Xcount
        Hc = Hc + lambda(i)*H{i};
    end
    Hc = Hc/sum(lambda);
    
    % Object function value and the relative error
    obj = 0;
    for i = 1:Xcount
        obj_body = norm(X{i} - H{i}*H{i}', 'fro')^2;
        obj_consensus = norm(H{i} - Hc, 'fro')^2;
        obj = obj + obj_body + lambda(i)*obj_consensus;
   
%         errX = mean(mean(abs(obj_body)))/mean(mean(X{i}));
%         errH = mean(mean(abs(obj_consensus)))/mean(mean(H{i}));
%         err = errX + errH;
    end
    
    
%     if rem(iter,5) == 0
%         fprintf([sprintf('Iter = '),int2str(iter),...
%             sprintf('\t Object function value = '),num2str(obj),...
%             sprintf('\t relative error = '),num2str(err), '\n']);
%     end
    
    Delta = obj_old - obj;
    obj_old = obj;
    objValue = [objValue, obj];
    
    if Delta < 1e-6
        break;
    end
end

end