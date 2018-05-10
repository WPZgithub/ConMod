function [ HPI ] = setSimilarity( modules )
moduleCounts = length(modules);
HPI = zeros(moduleCounts);
for i = 1:(moduleCounts-1)
    module_i = modules{i};
    for j = (i+1):moduleCounts
        module_j = modules{j};
        if ((~isempty(module_i)) && (~isempty(module_j)))
            HPI(i, j) = length(intersect(module_i,module_j))/min(length(module_i),length(module_j));
        end
    end
end
end