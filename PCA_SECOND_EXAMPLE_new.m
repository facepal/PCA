%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 2: Principal Component Analysis (PCA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear the command window and workspace
clc
clear all

% Set the desired explained variance (e.g., 0.95 for 95%)
explained_variance = 0.95;

% DATA SET IS THE MATRIX REPRESENTING THE IMAGE 'flowers.jpg'

% Read the image 'flowers.jpg'
A = imread('flowers', 'jpg');

% Display the original image in the first subplot
subplot(2,2,1)
imshow(A), axis off;
size(A); % Displays the size of A: NX x NY x 3 (color image with 3 RGB values per pixel)
title('Original Color Image')

% Convert the image to grayscale and cast to double precision
A = double(rgb2gray(A));
title('Original')

% Display the grayscale image in the third subplot
subplot(2,2,3)
imagesc(A), axis off; % Display the grayscale image
colormap gray; % Set colormap to gray
title('Grayscale Image')

% Assign the grayscale image to X for further processing
X = A;
aux = size(X);

N = aux(1); % Number of rows (pixels in vertical dimension)
m = aux(2); % Number of columns (pixels in horizontal dimension)

% Compute the mean of each column (mean over the rows)
media = mean(X)'; 

% Center the data by subtracting the mean from each column
Xc = zeros(size(X)); % Initialize the centered data matrix
for i = 1:m
    Xc(:, i) = X(:, i) - media(i);
end

% Compute the covariance matrix of the centered data
M1 = Xc' * Xc;

% Compute eigenvalues (DC) and eigenvectors (VPC) of the covariance matrix
% Compute up to 666 largest eigenvalues and corresponding eigenvectors
[VPC, DC] = eigs(M1, min(666, m)); 

% Determine the number of principal components required to achieve the desired explained variance
k = 0;
var_expl = 0;

% Sum the eigenvalues until the cumulative explained variance reaches the threshold
total_variance = trace(DC);
for i = 1:size(DC, 1)
    var_expl = var_expl + DC(i, i) / total_variance;
    k = k + 1;
    if var_expl >= explained_variance
        break;
    end
end

% Zero out the eigenvectors beyond the required number k
VPC(:, k+1:end) = 0;

% Retain only the necessary eigenvectors
VPC_util = VPC;

% Display the number of principal components needed and the rank of the covariance matrix
disp('Number of principal components to achieve explained variance:')
disp(k)
disp('Rank of covariance matrix M1:')
disp(rank(M1))

%pause % Pause execution to allow inspection of variables if needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct the data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Projection of data onto principal directions')

% Project the centered data onto the principal components
% Each row of pca_X contains the projection of a sample onto the principal components
pca_X = Xc * VPC_util;

% Display the size of the projected data matrix
size(pca_X);

% Reconstruct the data from the principal components
% Multiply the projected data by the transpose of the eigenvectors to map back to original space
rec_X = pca_X * VPC_util';

% Add the mean back to the reconstructed data to obtain the final reconstructed image
XX = rec_X + repmat(media', N, 1);

% Display the reconstructed image in the fourth subplot
subplot(2,2,4)
imagesc(XX), axis off;
colormap gray; % Ensure the colormap is set to gray
title('Reconstructed Image')
