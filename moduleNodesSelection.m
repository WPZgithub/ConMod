function [ modulesFinal ] = moduleNodesSelection( Hc, xita )
%  A soft node selection procedure from the consensus factors to assign the module members 
% and then truing the  modules to obtain more accurate results
%
% INPUT:
%   Hc: the consensus factor matrix
%   xita: the parameter for selecting nodes
%
% OUTPUT:
%   modulesFinal: a cell contains the final result modules
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

[N, K] = size(Hc);
candidateModules = cell(K, 1);
moduleSignal = zeros(K, 1);
H_mean = mean(Hc);
H_std = std(Hc, 0, 1);
for k = 1:K
    candidateModules{k} = find(Hc(:, k) > H_mean(k) + xita*H_std(k)); % Z-score>=t
    moduleSignal(k) = mean(Hc(candidateModules{k}, k));
end

HPI = setSimilarity( candidateModules );
modulesFinal = candidateModules;

for i = 1:size(HPI, 1)-1
    for j = (i+1):size(HPI, 2)
        if HPI(i,j)>0.5
            [Y, I] = max([moduleSignal(i), moduleSignal(j)]);
            if I == 1
                modulesFinal{j} = [];
                moduleSignal(j) = 0;
                HPI(j, :) = zeros(1, size(HPI, 2));
                HPI(:, j) = zeros(size(HPI, 1), 1);
            else
                modulesFinal{i} = [];
                moduleSignal(i) = 0;
                HPI(i, :) = zeros(1, size(HPI, 2));
                HPI(:, i) = zeros(size(HPI, 1), 1);
            end
        end
    end
end

i = 1;
while i ~= length(modulesFinal)+1
    if isempty(modulesFinal{i}) || (length(modulesFinal{i})<5)
        modulesFinal(i) = [];
        i = i - 1;
    end
    i = i + 1;
end

end