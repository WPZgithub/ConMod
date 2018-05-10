===============================================
    ConMod: Identifying conserved functional modules in multiple networks
===============================================

#----------------------The main function---------------------------
ConMod.m 	% The implement of ConMod method
>>  [ modulesfinal ] = ConMod( multiNetworks, N, K, lambda, xita, maxIter )
% INPUT:
%       multiNetworks: a cell contains multiple networks, each is  presented by a sparse matrix or a full matrix with N nodes
%       N: the number of all nodes
%       K: the number of hidden factors
%       lambda: a vector containing the parameters for balancing the relative
%		       weight among different views
%       xita:  the parameter for selecting nodes
%       maxIter:  the maximum number of iterations for multi-view NMF
%
% OUTPUT:
%       modulesfinal: a cell contains the final conserved modules


---------------The code for generating the synthetic datasets---------
syn_dataset_common.m		% Conserved modules have the same size and are common to a given set of networks
syn_dataset_overlap.m		% Conserved modules are present only in a subset of networks and they are the overlapping parts of specific modules across different networks

--------------The code for evaluation measures-------------
evaluation.m		% Compute the performance measures (TPR, FPR, Accuracy and MCC)

------------------------------------------------------------------

If any questions, please contact wangpeizhuo_37@163.com.