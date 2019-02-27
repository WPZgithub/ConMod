function [ Strength, Participation ] = featureNets( multiNetworks, N )
% Compute two feature metrices from multiple networks
%
% INPUT:
%   multiNetworks : a cell contains multiple networks, each of which is
%                   presented by edgelist format or a full matrix
%                   with N nodes
%   N : the number of all nodes
%
% OUTPUT:
%   Strength : N x N matrix for Connection Strength
%   Participation : N x N matrix for Participation Coefficient
%
% Peizhuo Wang (wangpeizhuo_37@163.com)

network_count = length(multiNetworks);
[n, m] = size(multiNetworks{1});

Strength = zeros(N);
temp = zeros(N);
A = zeros(N);
if (m <= 3) % Edgelist format
    for k = 1:network_count
        theMatrix = multiNetworks{k};
        [edge_count, col_count] = size(theMatrix);
        weight_max = 1;
        weight_min = 0;
        if (col_count == 3)
            weight_max = max(multiNetworks{k}(:,3));  
        end
        for e = 1:edge_count
            ii = theMatrix(e, 1);
            jj = theMatrix(e, 2);
            if (col_count == 3) % weighted network
                weight = abs(theMatrix(e, 3));
            else
                weight = 1; % unweighted network
            end
           
            Strength(ii, jj) = Strength(ii ,jj) + weight;
            Strength(jj, ii) = Strength(ii, jj);
            
            weight_1 = (weight - weight_min) / (weight_max-weight_min);
            weight_2 = 1/(1+exp(log(9999)-2*log(9999)*weight_1));
            if (weight_1 <= 0.3)
                weight_2 = 0;
            end
            
            A(ii, jj) = A(ii, jj) + weight_2;
            A(jj, ii) = A(ii, jj);
            temp(ii, jj) = temp(ii, jj) + weight_2^2;
            temp(jj, ii) = temp(ii, jj);
        end
    end
elseif (m == n) % Full matrix format
    N = n;
    for k = 1:network_count        
        % The edge weight is transformed using a logistic function, such that for the
        % element less than 0.3, we make it close to 0; for the element more than
        % 0.6, we make it close to 1.
        weight_max = max(max(multiNetworks{k}));
        weight_min = 0;
        matrix_weight_1 = (multiNetworks{k} - ones(N)*weight_min) ./ (weight_max-weight_min);
        matrix_weight_2 = 1./(1+exp(log(9999)-2*log(9999)*matrix_weight_1));

        matrix_weight_2(matrix_weight_1 <= 0.3) = 0;
        A = A + matrix_weight_2;
        temp = temp + matrix_weight_2.^2;
        Strength = Strength + multiNetworks{k};
    end
end
    
Participation = (network_count/(network_count-1)) * (1-(temp./(A.^2)));
Participation(isinf(Participation)) = 0;
Participation(isnan(Participation)) = 0;
Participation = Participation - diag(diag(Participation)); % The diagonal is 0

Strength = A./network_count;
Strength = Strength - diag(diag(Strength)); % The diagonal is 0

end