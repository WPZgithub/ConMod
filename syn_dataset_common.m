function [ dataset, realLabels ] = syn_dataset_common(alpha, isSaveToFiles, path)
% Synthetic dataset #1
% Conserved modules have the same size and are common to a given set of
% networks.
%
% INPUT:
%   alpha: The probability of the edge connected inside a module
%   isSaveToFiles: a flag for deciding whether to store the results to
%                  files (only for sparse matrix format)
%   path: output file path (it must be denoted if 'isSaveToFiles' is true)
%
% OUTPUT:
%   dataset: the generated dataset
%   realLabels: labels indicating the conserved module to which each point is allocated.
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

N = 500; % Number of nodes
M = 30; % Number of networks
dataset = cell(M, 1);
C_size = 80; % Size of each cluster
C1 = (1:C_size);
C2 = C1 + C_size;
C3 = C2 + C_size;
C4 = C3 + C_size;
C5 = C4 + C_size;
C = {C1, C2, C3, C4, C5};
realLabels = C;
M_C = length(C) * 5; % Number of networks with defined patterns

p_in = alpha;
p_out = 0.05;
if (p_out*(N-C_size))>(p_in*(C_size-1))
    p_out = (p_in*(C_size-1)) / (N-C_size); % (N-n)*p_out < (n-1)*p_in
end

mc = 1;
C_temp = {};
for m = 1:M
    % % Background network
    W = unifrnd(0,1,N,N);
    WW = zeros(N);
    WW(W < p_out) = 1;
    
    if (m <= M_C)
        if (mod(m-1, 5) == 0)
            C_temp = [C_temp, C{mc}];
            mc = mc + 1;
        end
        for i = 1:length(C_temp)
            % inside the cluster
            C_in = C_temp{i};
            WW1 = zeros(C_size);
            WW1(W(C_in, C_in) < p_in) = 1; 
            WW(C_in, C_in) = WW1;
        end
    else
        % One random module in each of the last 5 networks
        s = randperm(N);
        C_in = s(1:80);
        WW1 = zeros(length(C_in));
        WW1(W(C_in, C_in) < p_in) = 1; 
        WW(C_in, C_in) = WW1;
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
    dataset{m} = WW;
    
    % Save to file. Sparse matrix. Delete the edges with weight less than 0.3.
    if isSaveToFiles
        fp = fopen([path, 'network_', num2str(m)],'wt');
        for i=1:(size(WW,1)-1)
            for j=(i+1):size(WW,1)
                if (WW(i, j) >= 0.25)
                    fprintf(fp, '%d\t%d\t%f\n', i, j, WW(i,j));
                end
            end
        end
        fclose(fp);      
        fp = fopen([path, 'labels.txt'],'wt');
        for i=1:length(realLabels)
            theLabel = realLabels{i};
            for j=1:length(theLabel)
                fprintf(fp, '%d\t', theLabel(j));
            end
            fprintf(fp, '\n');
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