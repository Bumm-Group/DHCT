function [xScaled, yScaled, fastOffset, slowOffset] = ImCoords(RectImage, solution_vector, XScale, Period, OrigYPixels, SlowScan, FastScan)

%     ver 9/16/2016
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
% ImCoords calculates new image coordinates for an image with the same size
% as RectImage using the parameters in solution_vector, as determined by
% molecule_local_distortion. These newly calculated coordinates can be
% applied to generate a resampled image with ImReshape.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         ImCoords.m (this script)
%            -drift_models.m
%
%     Dependencies
%         drift_models.m
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

% Condition input data

if ~iscell(RectImage)
    RectImage = {RectImage};
    num_images = 1;
else
    num_images = length(RectImage);
end

if SlowScan == 0 || SlowScan == 1
    dimSlow = size(RectImage{1}, 2); % Image size in pixels
    dimFast = size(RectImage{1}, 1) + 2; 
else
    dimSlow = size(RectImage{1}, 1); % Image size in pixels
    dimFast = size(RectImage{1}, 2) + 2;
end

% Define useful constants

nm_to_m = 1e-9;

scanspeed = XScale / Period / nm_to_m;
slowscan = scanspeed / OrigYPixels / 2;
dist_to_pix = nm_to_m / XScale;

fast_linetime = Period * OrigYPixels;
image_time = Period * 2 * OrigYPixels.^2;

% Define corrections as functions of time

syms t

[xDrift(t), yDrift(t), h(t)] = drift_models('Symfun', solution_vector);

% Integrate velocity to determine position offset as a function of time

inth = int(h, t);
intx = int(xDrift, t);
inty = int(yDrift, t);

for i = 1:num_images
    
    % Determine fast-scan time
        
    fast_point = (0:dimFast - 1) / dist_to_pix;
    
    time_h = fast_point / scanspeed;
    if FastScan(i) == 1 || FastScan(i) == 2
        time_h = fliplr(time_h);
    end
    
    % Find the value of inth, intx and inty at t = 0
    
    if FastScan(i) == 1 || FastScan(i) == 2
        h0 = double(inth(fast_linetime));
    else
        h0 = double(inth(0));
    end
    
    if SlowScan == 1 || SlowScan == 2
        x0 = double(intx(image_time));
        y0 = double(inty(image_time));
    else
        y0 = double(inty(0));
        x0 = double(intx(0));
    end
      
    % Determine slow scan time
    
    slow_point = (1:dimSlow) / dist_to_pix;
    
    time_d = slow_point / slowscan;
    
    if SlowScan == 1 || SlowScan == 2
        time_d = fliplr(time_d);
    end
    
    % Calculate how much each pixel is offset
    
    Xpos = (double(intx(time_d)) - x0) * dist_to_pix;
    Ypos = (double(inty(time_d)) - y0) * dist_to_pix;
    Hpos = (double(inth(time_h)) - h0) / Period;
    
    % Choose appropriate offset direction
    
    if FastScan(i) == 1 || FastScan(i) == 2
        Hpos = -Hpos;
    end
    
    if SlowScan == 1 || SlowScan == 2
        Xpos = -Xpos;
        Ypos = -Ypos;
    end
    
    % Generate old x and y pixel coordinates
    
    [xgrid, ygrid] = meshgrid(fast_point * dist_to_pix, slow_point * dist_to_pix);
    
    % Generate x and y pixel coordinate offsets
    
    xScaled{i} = (xgrid - (ones(dimFast, 1) * Xpos).' - ones(dimSlow, 1) * Hpos); %#ok<*AGROW>
    yScaled{i} = (ygrid - (ones(dimFast, 1) * Ypos).');
    
    fastOffset(i) = -min(xScaled{i}(:)) + 1;
    slowOffset(i) = -min(yScaled{i}(:)) + 1;
    
    % Add offset to old pixel locations to get new pixel locations
    
    xScaled{i} = xScaled{i} + fastOffset(i);
    yScaled{i} = yScaled{i} + slowOffset(i);
    
end