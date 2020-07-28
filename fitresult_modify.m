function fitresult_out = fitresult_modify(fitresult, solution_vector, XScale, Period, OrigYPixels, SlowScan, FastScan, xOffset, yOffset) 

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
% fitresult_modify modifies the coordinates of image features in the
% fitresult data structure in the same way that the image pixel coordinates
% are modified in ImCoords. It uses the drift solution determined by
% molecule_local_distortion.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         fitresult_modify.m (this script)
%            -drift_models.m
%
%     Dependencies
%         drift_models.m
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

fitresult_out = fitresult;

% Define useful constants

tPts = fitresult.Time;
tPts_h = fitresult.Fast_Time;
nm_to_m = 1e-9;

fast_linetime = Period * OrigYPixels;
image_time = Period * 2 * OrigYPixels.^2;

dist_to_pix = nm_to_m / XScale;

% Define corrections as functions of time

syms t

[xDrift(t), yDrift(t), h(t)] = drift_models('Symfun', solution_vector);

% Integrate velocity to determine position offset as a function of time

inth = int(h, t);
intx = int(xDrift, t);
inty = int(yDrift, t);

% Find the value of inth, intx and inty at t = 0

if FastScan == 1 || FastScan == 2
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

% Calculate how much each feature is offset

for i = 1:length(tPts)
    h_t(i) = (double(inth(tPts_h(i))) - h0) / Period ;
    
    x_t(i) = (double(intx(tPts(i))) - x0) * dist_to_pix;
    y_t(i) = (double(inty(tPts(i))) - y0) * dist_to_pix;
end

% Choose appropriate offset direction

if FastScan == 1 || FastScan == 2
    h_t = -h_t;
end

if SlowScan == 1 || SlowScan == 2
    x_t = -x_t;
    y_t = -y_t;
end

% Add offset to old feature locations to get new feature locations

if SlowScan == 2 || SlowScan == 3
    fitresult_out.X = (fitresult.X - h_t - x_t) + xOffset;
    fitresult_out.Y = (fitresult.Y - y_t) + yOffset;
else
    fitresult_out.Y = (fitresult.Y - h_t - x_t) + xOffset;
    fitresult_out.X = (fitresult.X - y_t) + yOffset;
end
drawnow;