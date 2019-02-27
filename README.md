# ConMod

This algorithm is used for identifying conserved functional modules in multiple networks, as described in:
> Feature related multi-view nonnegative matrix factorization for identifying conserved functional modules in multiple biological networks. Peizhuo Wang, Lin Gao, Yuxuan Hu and Feng Li. BMC Bioinformatics, 2018, 19(1): 394.

This algorithm is implemented primarily in Matlab 2015b.

## Usage

### Input format

The code takes a series of networks as an input. These networks must be stored in a variable of type *cell* in matlab. Each network can be represented by the following two types:

1. adjacency matrix
2. edge list, e.g:
```
1 2 0.62
1 3 0.88
2 9 0.14
...
```

These codes are for generating synthetic datasets:

- `syn_dataset_common.m`: Conserved modules have the same size and are common to a given set of networks.
- `syn_dataset_overlap.m`: Conserved modules are present only in a subset of networks and they are the overlapping parts of specific modules across different networks.

### Main functions

- `ConMod.m`: The implementation of the ConMod algorithm.
```
function modulesfinal = ConMod(multiNetworks, N, K, lambda, xita, maxIter)
% INPUT:
%   multiNetworks: a cell contains multiple networks, each of which is presented by edgelist format or a full matrix with N nodes
%   N: the number of all nodes
%   K: the number of hidden factors
%   lambda: a vector which contains the parameters for balancing the relative weight among different views
%   xita:  the parameter for selecting nodes
%   maxIter:  the maximum number of iterations for multi-view NMF
%
% OUTPUT:
%   modulesfinal: a cell which contains the final conserved modules
```

- `featureNets.m`: Compute two feature matrices which characterize the multiple networks.
```
function [Strength, Participation] = featureNets(multiNetworks, N)
% INPUT:
%   multiNetworks: a cell contains multiple networks, each is  presented by a sparse matrix or a full matrix with N nodes
%   N: the number of all nodes
%
% OUTPUT:
%   Strength : N x N matrix for Connection Strength
%   Participation : N x N matrix for Participation Coefficient
```

- `multiViewNMF.m`: Multi-view non-negative symmetric matrix factorization.
```
function [H, Hc, objValue] = multiViewNMF(X, K, lambda, maxIter)
% INPUT:
%   X: a cell which contains symmetric matrices
%   K: the number of hidden factors
%   lambda: a vector which contains the parameters for balancing the relative weight among different views
%   maxiter: the maximum number of iterations
%
% OUTPUT:
%   H: a cell containing factor matrices for all views
%   Hc: the result consensus factor matrix
%   objValue: the value of objective function
```

- `moduleNodesSelection.m`: Assigning the module members by a soft node selection procedure and then truing the modules to obtain more accurate results
```
function modulesFinal = moduleNodesSelection(Hc, xita)
% INPUT: 
%   Hc: the consensus factor matrix
%   xita: the parameter for selecting nodes
%
% OUTPUT:
%   modulesFinal: a cell which contains the final result modules
```

If you have any questions, please contact `wangpeizhuo_37@163.com`.