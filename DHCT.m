function out = DHCT(XScale, Period, SlowScanDir, FastScanDir, spacing, varargin)

%     ver 9/16/2016
%
%     copyright (c) 2016 Mitchell P. Yothers & Lloyd A. Bumm
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
% DHCT is a script that runs all of the stripts required for STM
% image correction due to hysteresis, creep and thermal drift in order,
% returning an output structure with a variety of useful data including the
% corrected image and the locations of the image features. It requires as
% inputs the space between adjacent pixels in meters (XScale), the time
% between subsequent pixels in seconds (Period), Strings for the slow and
% fast scan directions ('Up', 'Down', 'Left', or 'Right'), and the
% nearest-neighbor spacing in meters for features on your surface determined
% by some other technique, like XRD or LEED, and the image that you want to
% correct. A pair of images can be corrected by including both images and
% their fast scan directions in the arguments - FastScanDir should then be
% a cell array with both direction strings in it. Logical region masks the
% same size as the input images can additionally be supplied, only
% positively-masked areas will be used for the image correction. 
% 
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m (this script)
%         data_mask_select.m
%            -find_nearest_neighbors.m
%            -position_analysis.m
%         feat_position.m
%         fitresult_modify.m
%            -drift_models.m
%         find_molecules.m
%            -fmgaussfit_reduced.m
%         ImCoords.m
%            -drift_models.m
%         ImReshape.m
%         molecule_local_distortion.m
%            -drift_models.m
%         time_assign.m
%
%     Dependencies
%         data_mask_select.m
%         find_nearest_neighbors.m
%         position_analysis.m
%         feat_position.m
%         fitresult_modify.m
%         drift_models.m
%         find_molecules.m
%         fmgaussfit_reduced.m
%         ImCoords.m
%         ImReshape.m
%         molecule_local_distortion.m
%         time_assign.m
% 
%     Dependents (*called by) 
%         none 
% 

% If show_images is set to 1, Matlab will generate figures showing some
% useful data during its progress on correcting the image. If set to 0,
% these will be suppressed. 

show_images = 1;

% Condition input data

Image{1} = varargin{1};

if nargin == 7
    if islogical(varargin{2})
        Mask{1} = varargin{2}; %#ok<*NASGU>
        num_images = 1;
    else
        Image{2} = varargin{2};
        Mask{1} = logical(Image{1} .* 0 + 1);
        Mask{2} = logical(Image{2} .* 0 + 1);
        num_images = 2;
    end
elseif nargin > 7
    Image{2} = varargin{2};
    Mask{1} = varargin{3};
    Mask{2} = varargin{4};
    num_images = 2;
else
    num_images = 1;
end

% Turn directional strings into a number for more convenient use in
% subsequent code

if strcmpi(SlowScanDir, 'Right');
    SlowScan = 0;
elseif strcmpi(SlowScanDir, 'Left');
    SlowScan = 1;
elseif strcmpi(SlowScanDir, 'Up');
    SlowScan = 2;
elseif strcmpi(SlowScanDir, 'Down');
    SlowScan = 3;
else
    error('Slow Scan Direction not recognized.');
end

for i = 1:num_images
    if strcmpi(FastScanDir{i}, 'Right');
        FastScan(i) = 0; %#ok<*AGROW>
    elseif strcmpi(FastScanDir{i}, 'Left');
        FastScan(i) = 1;
    elseif strcmpi(FastScanDir{i}, 'Up');
        FastScan(i) = 2;
    elseif strcmpi(FastScanDir{i}, 'Down');
        FastScan(i) = 3;
    else
        error('Fast Scan Direction not recognized.');
    end
end

% Number of pixels in the slow scan direction on the original image is
% important for a variety of things during the correction

if SlowScan == 2 || SlowScan == 3
    out.SlowPixels = size(Image{1}, 2);
else
    out.SlowPixels = size(Image{1}, 1);
end

% Turn the nearest neighbor spacing in meters into pixels

spacing_guess = spacing / XScale;

% This code prints out what it is doing in the main Matlab window between 
% each large chunk of the process. For more information on each of these 
% pieces, refer to their comments.

fprintf('\nFinding feature locations...\n')

for i = 1:num_images
    fit_data{i} = find_molecules(Image{i}, spacing_guess, show_images);
end

fprintf('\nAssigning feature times...\n')

out.fit_data = time_assign(fit_data, XScale, Period, out.SlowPixels, SlowScan, FastScan);

fprintf('\nFinding nearest neighbors...\n')

for i = 1:num_images
    out.mask_fit_data{i} = data_mask_select(out.fit_data{i}, Mask{i}, XScale, spacing, show_images);
end

% Show feature locations on each image, as selected by the input masks.

if show_images
    for i = 1:num_images
        figure;
        imshow(Image{i}, []);
        hold on;
        scatter(out.mask_fit_data{i}.X, out.mask_fit_data{i}.Y, '.', 'red');
    end
end

fprintf('\nFinding Distortion from nearest neighbors...\n')
 
out.feat_pos_nn = feat_position(out.mask_fit_data, XScale, SlowScan);

[out.solution_vector, out.driftX, out.driftY, out.scaleFactor] ...
    = molecule_local_distortion(out.mask_fit_data, XScale, Period, ...
    out.SlowPixels, spacing, out.feat_pos_nn, 'Both', show_images);

fprintf('\nCorrecting Image...\n')

[out.xCoords, out.yCoords, out.xOffset, out.yOffset] = ImCoords(Image, out.solution_vector, XScale, Period, out.SlowPixels, SlowScan, FastScan);

for i = 1:num_images
    out.ReshapedImage{i} = ImReshape(Image{i}, SlowScan, out.xCoords{i}, out.yCoords{i}, show_images);
end

fprintf('\nCorrecting Molecular Indices...\n')

for i = 1:num_images
    out.mask_drift_data{i} = fitresult_modify(out.mask_fit_data{i}, out.solution_vector, ...
        XScale, Period, out.SlowPixels, SlowScan, FastScan(i), out.xOffset(i), out.yOffset(i));
end

fprintf('\nDone!\n');