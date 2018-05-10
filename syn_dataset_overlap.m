function [dataset, realLabels, lables_specific] = syn_dataset_overlap(alpha, isSaveToFiles, path)
% Synthetic dataset #2
% Conserved modules are present only in a subset of networks and they are
% the overlapping parts of specific modules across different networks.
%
% INPUT: 
%   alpha: the probability of the edge connected inside a module
%   isSaveToFiles: a flag for deciding whether to store the results to
%                  files (only for sparse matrix format)
%   path: output file path (it must be denoted if 'isSaveToFiles' is true)
%
% OUTPUT:
%   dataset: the generated dataset
%   realLabels: labels indicating the conserved module to which each point is allocated.
%   lables_specific:lables of each specific module on each network
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

N = 500; % Number of nodes
M = 15; % Number of networks
dataset = cell(M, 1);
C_common_1 = 1:50; % Common part 1
C_rest_1 = setdiff((1:N), C_common_1);
C_common_2 = 401:440; % Common part 2
C_rest_2 = setdiff((1:N), C_common_2);
C_rest = intersect(C_rest_1, C_rest_2);
realLabels = {C_common_1, C_common_2};
lables_specific = zeros(N, M);

p_in = alpha;
k = 0;
for m = 1:M
    if (m <= 5)
        s = C_rest_1(randperm(length(C_rest_1)));
        increment = 10 * (11-m);
        C_specific = {[C_common_1, s(1:increment)]}; % 150:10:110
    elseif ((5 < m) && (m <= 10))
        s_1 = C_rest(randperm(length(C_rest)));
        C_specific_1 = [C_common_1, s_1(1:(10 * (11-m)))]; % 100:10:60
        s_1_rest = setdiff(s_1, C_specific_1);
        s_2 = s_1_rest(randperm(length(s_1_rest)));
        C_specific_2 = [C_common_2, s_2(1:(5 * (m-2)))]; % 60:5:80
        C_specific = {C_specific_1, C_specific_2};  
    else
        s = C_rest_2(randperm(length(C_rest_2)));
        increment = 5 * (m-2);
        C_specific = {[C_common_2, s(1:increment)]}; % 85:5:105
    end
    
    p_out = 0.05;
    if (p_out*(N-60))>(p_in*(60-1))
        p_out = (p_in*(60-1)) / (N-60); % (N-n)*p_out < (n-1)*p_in
    end
    % Background network
    W = unifrnd(0,1,N,N);
    WW = zeros(N);
    WW(W < p_out) = 1;
    k = k + 1; 
    
    % inside the module
    for i = 1:length(C_specific)
        lables_specific(C_specific{i}, m) = i;
        WW1 = zeros(length(C_specific{i}));
        WW1(W(C_specific{i}, C_specific{i}) < p_in) = 1;
        WW(C_specific{i}, C_specific{i}) = WW1;
    end
    
    % Gaussian noise, sigma=0.1 or 0.15
    WW_tril = tril(WW);
    E = normrnd(0.25, 0.1, N, N);
    WW_0 = WW_tril + tril(E); % X0+E
    WW_0(WW_tril == 1) = 0;
    WW_0(WW_0 < 0) = 0;
    WW_0(WW_0 > 1) = 1;
    WW_1 = WW_tril - tril(E); % X1-E
    WW_1(WW_tril == 0) = 0;
    WW_1(WW_1 > 1) = 1;
    WW_1(WW_1 < 0) = 0;
    WW = tril(WW_1 + WW_0);

    WW = WW - diag(diag(WW));
    WW = tril(WW) + tril(WW)';
    dataset{k} = WW;

    % Save to file. Sparse matrix. Delete the edges with weight less than 0.3.
    if isSaveToFiles
        fp = fopen([path, '\network_', num2str(m), '.txt'],'wt');
        for i=1:(size(WW,1)-1)
            for j=(i+1):size(WW,1)
                if (WW(i, j) >= 0.25)
                    fprintf(fp, '%d\t%d\t%f\n', i, j, WW(i,j));
                end
            end
        end
        fclose(fp);
    end
end
if isSaveToFiles
    flist = fopen([path, 'networklist.txt'],'wt');
    for m = 1:M
        fprintf(flist, 'network_%d\n', m);    
    end
    fclose(flist);
end

end