function [ pvalues_modulePerNet, FDR ] = significantModules( modules, multiNetworks, N, permu_times )
% Using a permutation test to assess the significance of functional modules across multiple networks.
% This allows identifying the specific conditions where each module is detected.
% 
% INPUT:
%   modules:
%   multiNetworks: a cell containing adjacency matrices of multiple
%   networks
%   N: the total number of nodes in multiple networks
%   permu_times: the number of permutation times(default value is 1000)
%
% OUTPUT:
%   pvalues_modulePerNet: individual p-value of each module in each network
%   FDR: Benjamin-Hochberg adjusted p-values
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

%% Initialization
moduleCounts = length(modules);
networkCounts = length(multiNetworks);
nodeCounts = N;
[x,y] = size(multiNetworks{1});
if nargin == 3
    permu_times = 1000;
end

%% Permutation test
moduleDensity = zeros(moduleCounts, networkCounts);
clusterQuality = zeros(moduleCounts, networkCounts);
null_counts = zeros(moduleCounts, networkCounts, permu_times);
FDR = zeros(moduleCounts, networkCounts);
for m = 1:moduleCounts
    fprintf('Module %d/%d...\n', m, moduleCounts);
    theModule = modules{m};
    moduleN = length(theModule)*(length(theModule)-1)/2;
    
    for n = 1:networkCounts
        theNetwork = multiNetworks{n};
        
        % Compute the individual cluster quality and density
        if (y <= 3) % sparse matrix format
            isA = ismember(theNetwork(:, 1), theModule);
            isB = ismember(theNetwork(:, 2), theModule);
            indexEdges = isA & isB;
            if (y == 3) % weighted network
                inDen = sum(theNetwork(indexEdges, 3));
                moduleDensity(m, n) = inDen/moduleN;
                outDen = (sum(theNetwork(:, 3)) - inDen)/(nodeCounts-length(theModule));
                outDen(outDen==0) = 1e-6;
                clusterQuality(m, n) = moduleDensity(m, n)/outDen;
            else % unweighted network
                inDen = length(theNetwork(indexEdges));
                moduleDensity(m, n) = inDen/moduleN;
                outDen = (size(theNetwork, 1) - inDen)/(nodeCounts-length(theModule));
                outDen(outDen==0) = 1e-6;
                clusterQuality(m, n) = moduleDensity(m, n)/outDen;
            end
        elseif (x == y) % full matrix format
            inDen = sum(sum(theNetwork(theModule, theModule)))/2;
            outDen = (sum(sum(theNetwork))/2 - inDen)/(nodeCounts-length(theModule));
            outDen(outDen==0) = 1e-6;
            moduleDensity(m, n) = inDen/moduleN;
            clusterQuality(m, n) = moduleDensity(m, n)/outDen;
        end
        
        % Compute the random individual cluster quality and density
        % And compute individual p-value
        for t = 1:permu_times
            randnum = randperm(nodeCounts);
            randModule = randnum(1:length(theModule));
            
            if (y <= 3) % sparse matrix format
                isA = ismember(theNetwork(:, 1), randModule);
                isB = ismember(theNetwork(:, 2), randModule);
                indexEdges = isA & isB;
                if (y == 3) % weighted network
                    inDen = sum(theNetwork(indexEdges, 3));
                    randModuleDensity = inDen/moduleN;
                    outDen = (sum(theNetwork(:, 3)) - inDen)/(nodeCounts-length(theModule));
                    outDen(outDen==0) = 1e-6;
                    randClusterQuality = randModuleDensity/outDen;
                else % unweighted network
                    inDen = length(theNetwork(indexEdges));
                    randModuleDensity = inDen/moduleN;
                    outDen = (size(theNetwork, 1) - inDen)/(nodeCounts-length(theModule));
                    outDen(outDen==0) = 1e-6;
                    randClusterQuality = randModuleDensity/outDen;
                end
                if (randClusterQuality >= clusterQuality(m, n))
                    null_counts(m, n, t) = 1;
                end
            elseif (x == y) % full matrix format
                randModule_matrix = theNetwork(randModule, randModule);
                inDen = sum(sum(randModule_matrix))/2;
                outDen = (sum(sum(theNetwork))/2 - inDen)/(nodeCounts-length(theModule));
                outDen(outDen==0) = 1e-6;
                ranModuleDensity = inDen/moduleN;
                randClusterQuality = ranModuleDensity/outDen;
                if (randClusterQuality >= clusterQuality(m, n))
                    null_counts(m, n, t) = 1;
                end
            end
        end
    end
end
pvalues_modulePerNet = sum(null_counts, 3)/permu_times;
for j = 1:networkCounts
    FDR(:,j) = mafdr(pvalues_modulePerNet(:,j),'BHFDR', true);
end

end