function fitresult = find_molecules(raw_image, spacing_guess, varargin)

%     ver 9/1/2016
%
%     copyright (c) 2016 Mitchell P. Yothers & Lloyd A. Bumm
%
%     This file is part of DHCT.
%
%     DHCT is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     DHCT is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with DHCT.  If not, see <http://www.gnu.org/licenses/>.
%  
%     contact: 
%     L. A. Bumm
%     440 W Brooks St
%     Homer L. Dodge Department of Physics & Astronomy
%     The University of Oklahoma
%     Norman, OK 73019
%     bumm@ou.edu
%
% SUMMARY
% find_molecules determines feature locations in an STM image with a cross-
% correlation using a 2D Gaussian as a kernel. The kernel size is chosen
% based on the expected feature size. It returns the structure fitresult,
% which stores feature parameters on a per-feature basis. Each feature has
% eight parameters stored in the fitresult structure - the (x, y, z) 
% coordinates of the base of the best-fit Gaussian, its standard deviation 
% in 2 orthogonal directions, its amplitude, the angle between the major 
% axis and the X axis, and an index for keeping track of features across
% multiple data strutures.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         find_molecules.m (this script)
%            -fmgaussfit_reduced.m
%
%     Dependencies
%         fmgaussfit_reduced.m
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%
if nargin > 2
    show_images = varargin{1};
else
    show_images = 0;
end

% Value of cross correlation that selects valid surface regions

threshold = 0.1;

% Calculate optimal kernel gaussian sigma and size for cross-correlation

radius = round(spacing_guess / 2);
sigma = radius / 2;
fprintf('Gaussian template sigma = %3.1f pixels\n', sigma);
diameter = 2 * radius + 1;
template = fspecial('gaussian', diameter, sigma);

% Generate cross-correlation image and crop it to original image size

[ImageX, ImageY] = size(raw_image);

raw_xcorr = normxcorr2(template, raw_image);
xcorr_image = raw_xcorr(radius + (1:ImageX), radius + (1:ImageY));

% Use watershed and threshold to break xcorr into distinct regions

watershed_mask = double(logical(watershed(-xcorr_image)));

xcorr_mask = and(watershed_mask, xcorr_image > threshold);

% Create image with only selected pixels for analysis

image_peak = zeros(ImageX, ImageY);
x_sub = (radius + 1):(ImageX - radius);
y_sub = (radius + 1):(ImageY - radius);
image_peak(x_sub, y_sub) = raw_image(x_sub, y_sub) .* xcorr_mask(x_sub, y_sub);

% Display selected features (Optional)

if show_images
    figure;
    imagesc(image_peak);
    colormap gray;
    axis image;
    drawnow;
end

% Separate selected features into regions

regions = bwconncomp(image_peak);
area_pix = cellfun(@length, regions.PixelIdxList);

fprintf('total peaks found: %d areas range from %d to %d\n', regions.NumObjects, min(area_pix), max(area_pix));

% Discard regions too big or too small to be a feature

median_area = median(area_pix(area_pix > 10));

filt_width = sqrt(10);

area_lb = median_area / filt_width;
area_ub = median_area * filt_width;

in_range = and(area_pix > area_lb, area_pix < area_ub);
num_regions = sum(in_range);

fprintf('peaks in range %d to %d pixels: %d\n', floor(area_lb), floor(area_ub), num_regions);

% Display feature size histogram (Optional)

if show_images
    area_size = 1:ceil(area_ub + median_area);
    h_area = hist(area_pix, area_size);

    figure;
    plot(area_size, h_area);
    hold on;
    scatter([area_lb, area_ub], [max(h_area)/2, max(h_area)/2]);
    drawnow;
end

% Make the set of pixels that belong to each region easily accessible

loc_regions = find(in_range);

for i = 1:sum(in_range)
    pixel_list{i} = regions.PixelIdxList{loc_regions(i)};
end

for i = 1:num_regions
        
    % Obtain pixel x, y coordinates from pixel index 

    fx{i} = mod(pixel_list{i} - 1, ImageX) + 1;    
    fy{i} = floor((pixel_list{i} - 1) / ImageX) + 1;

    % Pixel z coordinates are from the STM image
    
    fz{i} = raw_image(sub2ind(size(raw_image), fx{i}, fy{i}));
end

% Initialize parallel pool

gcp;

fprintf('Progress:\n');
fprintf([repmat(' ', 1, 20) '25%% v' repmat(' ', 1, 20) '50%% v' ...
    repmat(' ', 1, 20) '75%% v' repmat(' ', 1, 19) '100%% v' '\n\n']);

% Find best-fit Gaussians to feature pixels 

parfor i = 1:num_regions
    
    out(i, :) = fmgaussfit_reduced(fy{i}, fx{i}, fz{i}, 0); %#ok<*AGROW>

    if floor(mod(i, num_regions / 100)) == 0
        fprintf('\b%%\n');
    end
end

% Assign reasonable names to output structure variables

fitresult.Amplitude = out(:, 1).';
fitresult.Angle = out(:, 2).';
fitresult.SigmaX = out(:, 3).';
fitresult.SigmaY = out(:, 4).';
fitresult.X = out(:, 5).';
fitresult.Y = out(:, 6).';
fitresult.Z = out(:, 7).';
fitresult.Index = 1:length(fitresult.X);