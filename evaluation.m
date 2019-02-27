function [ TPR, FPR, Accuracy, MCC ] = evaluation(C_preticted, C_reference, N)
% Compute the performance measures (TPR, FPR, Accuracy and MCC)
%
% INPUT: 
%       C_preticted: a cell containing the preticted cluster results
%       C_reference: a cell containing the ground truth
%       N: total number of items 
% OUTPUT: 
%       TPR: the True Positive Rate
%       FPR: the False Positive Rate
%       Accuracy:
%       MCC: the Matthews Correlation Coefficient
%       I: Confusion Matrix
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

%% Confusion matrix
K_preticted = length(C_preticted);
K_reference = length(C_reference);
Ncount_P = 0;
Ncount_R = 0;
c_num = zeros(K_preticted+1, K_reference+1);
C = cell(K_preticted, K_reference);
for i = 1:K_preticted
    Ncount_P = Ncount_P + length(C_preticted{i});
    theRef_nonoverlap = [];
    for j = 1:K_reference
        C{i, j} = intersect(C_preticted{i}, C_reference{j});
        c_num(i, j) = length(C{i, j});
        if (i == 1)
            Ncount_R = Ncount_R + length(C_reference{j});
        end
        theRef_nonoverlap = union(theRef_nonoverlap, C{i, j});
    end
    c_num(i, j+1) = length(C_preticted{i}) - length(theRef_nonoverlap); % Background noise nodes
end
for j = 1:K_reference
    thePre_nonoverlap = [];
    for i = 1:K_preticted
        thePre_nonoverlap = union(thePre_nonoverlap, C{i, j});
    end
    c_num(i+1, j) = length(C_reference{j}) - length(thePre_nonoverlap); % Lost reference nodes
end

%% TP, FP, FN, TN for nodes pairs
TP = sum(sum(c_num(1:K_preticted, 1:K_reference).*(c_num(1:K_preticted, 1:K_reference)-1)/2));
FP = 0;
for j = 1:K_reference
    tempC = c_num(1:K_preticted, (j+1):(K_reference+1));
    FP = FP + sum(c_num(1:K_preticted,j).*sum(tempC, 2));
end
FN = 0;
for i = 1:K_preticted
    tempC = c_num((i+1):(K_preticted+1), 1:K_reference);
    FN = FN + sum(c_num(i,1:K_reference).*sum(tempC, 1));
end
TN = N*(N-1)/2 - (TP+FN+FP);
I = [TP, FP; FN, TN];

%% TPR, FPR, MCC
TPR = TP/(TP+FN);
FPR = FP/(FP+TN);
Accuracy = (TP+TN)/(TP+FP+TN+FN);
MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));

end