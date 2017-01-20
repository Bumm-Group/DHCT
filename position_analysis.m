function [iLoc, jLoc, dist] = position_analysis(fitresult, XScale)

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
% Given the structure of feature locations from find_molecules, determine
% the distance between every pair of features. 
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         data_mask_select.m
%            -position_analysis.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         data_mask_select.m (*) 
%

num_features = length(fitresult.Index);

nm_to_m = 1e-9;

iLoc = [];
jLoc = [];

% Initialize parallel pool

gcp;

% Determine indices of pairs of features to compute distance

parfor i = 1:num_features - 1
    iLoc = [iLoc, i * ones(1, num_features - i)]; %#ok<*AGROW>
    jLoc = [jLoc, i + 1:num_features];
end


if gpuDeviceCount > 0

    % Get coordinates from indices and send to GPU. 

    gx1 = gpuArray(fitresult.X(iLoc));
    gx2 = gpuArray(fitresult.X(jLoc));
    gy1 = gpuArray(fitresult.Y(iLoc));
    gy2 = gpuArray(fitresult.Y(jLoc));

    % Calculate distance from coordinates on the GPU. 

    rsq = (gx2 - gx1) .^ 2 + (gy2 - gy1) .^ 2;
    dist = gather(sqrt(rsq) * XScale / nm_to_m);
    
else
    
    % Calculate distance from coordinates on the CPU when no GPU.
    
    rsq = (fitresult.X(jLoc) - fitresult.X(iLoc)) .^ 2 + (fitresult.Y(jLoc) - fitresult.Y(iLoc)) .^ 2;
    dist = sqrt(rsq) * XScale / nm_to_m;

end
