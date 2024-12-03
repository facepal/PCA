%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 1: Principal Component Analysis (PCA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear the workspace and command window, and set numerical display format
format long
clc
clear all

% The dataset has 4 samples (observations), each with 6 attributes (variables)
% Each row represents a sample
% Set the desired explained variance (1 corresponds to 100%)
var_explicada = 1;

% Define the data matrix X (after transpose, samples are rows and attributes are columns)
X = [1 2 3 4;
     11 5 3 2;
     1 18 24 1;
     17 12 5 4;
     27 30 66 88;
     2 1 3 4;]';

% Get the size of the data matrix
aux = size(X);

N = aux(1); % Number of samples (rows)
m = aux(2); % Number of attributes (columns)

% Compute the mean of each attribute (column-wise)
media = mean(X);

% Center the data by subtracting the mean from each attribute
Xc = zeros(size(X)); % Initialize centered data matrix
for i = 1:m
    Xc(:, i) = X(:, i) - media(i);
end

% Compute the covariance matrix of the centered data
M1 = Xc' * Xc;

% Compute eigenvalues (DC) and eigenvectors (VPC) of the covariance matrix
% Columns of VPC are eigenvectors, diagonal elements of DC are eigenvalues
[VPC, DC] = eigs(M1, min(100, m)); % Compute up to 100 largest eigenvalues (limited by m)

% Determine the number of principal components required to explain the desired variance
k = 0;
var_expl = 0;

for i = 1:N
    if var_expl <= var_explicada
        var_expl = var_expl + DC(i, i) / trace(DC);
        k = k + 1;
    end
end

% Zero out the eigenvectors beyond the required number k
for i = k+1:m
    for j = 1:length(M1)
        VPC(j, i) = 0.0;
    end
end

% Retain only the necessary eigenvectors
VPC_util = VPC;

%pause % Pause execution to allow inspection of variables if needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct the data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Projection of data onto principal directions');

% Project the centered data onto the principal components
% Each row of pca_X contains the projection of a sample onto the principal components
pca_X = Xc * VPC_util;

% Compute projections of individual samples (optional)
amostra1_proj = Xc(1, :) * VPC_util; % Projection of sample 1
amostra2_proj = Xc(2, :) * VPC_util; % Projection of sample 2
amostra3_proj = Xc(3, :) * VPC_util; % Projection of sample 3
amostra4_proj = Xc(4, :) * VPC_util; % Projection of sample 4

% Display the size of the projected data matrix
size(pca_X);

% Reconstruct the data from the principal components
% Multiply the projected data by the transpose of the eigenvectors to map back to original space
rec_X = pca_X * VPC_util'; % Reconstructed centered data

% Add the mean back to the reconstructed data to obtain the final reconstructed data
XX = zeros(size(X)); % Initialize reconstructed data matrix
for i = 1:m
    for j = 1:N
        XX(j, i) = rec_X(j, i) + media(i);
    end
end

disp('Reconstructed data');
XX

disp('Original data');
X
