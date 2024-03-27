function Anomaly_Img_to_Mask_Conversion(filepath)
% Anomaly_Img_to_Mask_Conversion.m
% Function to convert an image of a material property distribution to a 
% binary mask
% 
% Input argument:
% filepath (string): path to image of material property
%
% Author: Matthieu Ruthven (matthieu.ruthven@uni.lu)
% Last modified: 27th March 2024

% Read file
A = imread(filepath);

% Convert RGB to greyscale
I = rgb2gray(A);

% Unique pixel values
unique_pix = unique(I);

% Preallocate vector
pix_count = zeros(size(unique_pix));

% Preallocate counter
counter = 1;

% For each value
for pix_val = unique_pix'

    % Calculate number of pixels
    pix_count(counter) = numel(find(I == pix_val));

    % Update counter
    counter = counter + 1;

end

% Sort pix_count in descending order
[B, idx] = sort(pix_count, 'descend');

% Pixel value of background
bkg = unique_pix(idx(1));

% Pixel value of majority of battery
btt = unique_pix(idx(2));

% Identify other regions of image
other = ((I ~= bkg) + (I ~= btt)) == 2;

% Identify battery boundary
BW = imfill(other, 'holes');
BW = bwperim(BW);

% Identify anomaly
anom = other - BW;
anom = imfill(anom, 'holes');
anom = bwareafilt(anom == 1, 1);

% Create PNG file
[filepath, name, ext] = fileparts(filepath); 
imwrite(anom * 1, fullfile(filepath, [name(1:(end-1)) 'Anomaly' ext]))