function [ modulesfinal ] = ConMod( multiNetworks, N, K, lambda, xita, maxIter )
% The implement of ConMod method for identifying conserved functional
% modules in multiple networks
%
% INPUT:
%   multiNetworks : a cell contains multiple networks, each of which is
%                   presented by edgelist format or a full matrix
%                   with N nodes
%   N: the number of all nodes
%   K: the number of hidden factors
%   lambda: a vector which contains the parameters for balancing the relative
%		   weight among different views
%   xita:  the parameter for selecting nodes
%   maxIter:  the maximum number of iterations for multi-view NMF
%
% OUTPUT:
%   modulesfinal: a cell which contains the final conserved modules
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

%% Calculting the feature matrices
disp('Calculating the strengh matrix and the uniformity matrix...')
[Strength, Distribution] = featureNets(multiNetworks, N);

%% Obtaining the candidate modules by multi-view NMF
disp('Obtaining candidate modules by multi-view NMF...')
disp(['K=', num2str(K)])
X = {Strength, Distribution};
[ H, Hc, objValue ] = multiViewNMF( X, K, lambda, maxIter );

%% Selecting nodes from the consensus factors
modulesfinal = moduleNodesSelection( Hc, xita );

end

