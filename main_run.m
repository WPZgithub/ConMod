clear
clc
%% Load multi-network data sets
disp('Processing the networks...')

% I: Load data from files
% nets = importdata('./networklist.txt');
% T = length(nets);
% multiNetworks = cell(T, 1);
% for i = 1:T
%     disp(['network: ', num2str(i)])
%     multiNetworks{i} = load(['./', nets{i}]);
% end
% realLabels = importdata('./labels.txt');

% II: Automatically generated synthetic datasets
%[multiNetworks, realLabels] = syn_dataset_common(0.5, false); % common type
[multiNetworks, realLabels, lables_specific] = syn_dataset_overlap(0.3, false); % overlap type
num_Nodes = 500;

%% One-step finding conserved functional modules
tic
K = 5;
lambda = [0.01, 0.05];
xita = 2;
maxIter = 50;
modules = ConMod( multiNetworks, num_Nodes, K, lambda, xita, maxIter );
runtime = toc;
disp(['Done.    Running time: ', num2str(runtime), ' sec.'])

%% Step-by-step finding conserved functional modules
% Calculting the feature networks
% tic
% disp('Calculating the strengh matrix and the uniformity matrix...')
% [Strength, Participation] = featureNets(multiNetworks, num_Nodes);
% 
% % Obtaining the candidate modules by multi-view NMF
% disp('Obtaining candidate modules by multi-view NMF...')
% K = 5;
% disp(['K=', num2str(K)])
% X = {Strength, Participation};
% lambda = [0.01, 0.05];
% [ H, Hc, objValue ] = multiViewNMF( X, K, lambda, 50 );
% 
% % Selecting nodes from the consensus factors
% xita = 1.5;
% modules_final = moduleNodesSelection( Hc, xita );
% runtime = toc;
% disp(['Running time: ', num2str(runtime), ' sec.'])

%% Module validation
% disp('Validation...')
% [ pvalues_modulePerNet, FDR2 ] = significantModules(modules_Merged, multiNetworks, num_Nodes);
% disp('Done.')

%% Clustering performance
[ TPR, FPR, Accuracy, MCC] = evaluation(modules, realLabels, num_Nodes);