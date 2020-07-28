function Image = ImReshape(RectImage, SlowScan, xScaled, yScaled, varargin)

%     ver 7/20/2020
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
% ImReshape resamples an input image for display from an unevenly spaced
% set of pixel coordinates, like those generated from ImCoords. 
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         ImReshape.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

% Set lowest image pixel value to zero

if min(RectImage(:)) ~= max(RectImage(:))
    RectImage = RectImage - min(RectImage(:));
end

% Add a row of zero pixels to the edges of the image to smoothly transition
% into the area with no data

if SlowScan == 0 || SlowScan == 1
    dimSlow = size(RectImage, 2); % Image size in pixels
    v = zeros(1, dimSlow);
    
    RectImage = cat(1, v, RectImage, v).'; 
    % Flip image orientation to match correction orientation
else
    dimSlow = size(RectImage, 1); % Image size in pixels
    v = zeros(dimSlow, 1);
    
    RectImage = cat(2, v, RectImage, v);
end

% Generate resampled image

[xPixels, yPixels] = meshgrid(1:floor(max(xScaled(:))), 1:floor(max(yScaled(:))));

Image = griddata(xScaled(:), yScaled(:), RectImage(:), xPixels, yPixels);

% Return image to correct orientation

if SlowScan == 0 || SlowScan == 1
    Image = Image.';
end

% Replace pixels which have no data with zeroes

Image(isnan(Image)) = 0;

% Show image (Optional)

if nargin > 4
    show_images = varargin{1};
    if show_images
        figure;
        imagesc(Image);
        axis image;
        colormap(gray);
    end
end