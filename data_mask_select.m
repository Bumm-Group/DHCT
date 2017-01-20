function data_out = data_mask_select(data, Mask, XScale, spacing, varargin)

%     ver 9/5/2016
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
% data_mask_select keeps only features from the data structure that are 
% in the logical mask supplied in Mask, and determines their nearest
% neighbors from the radial distribtuion function of the data. 
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         data_mask_select.m (this script)
%            -find_nearest_neighbors.m
%            -position_analysis.m
%
%     Dependencies
%         find_nearest_neighbors.m
%         position_analysis.m
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

% Define useful constants

if nargin > 4
    show_images = varargin{1};
else
    show_images = 0;
end

nm_to_m = 1e-9;

r_0 = spacing / nm_to_m;

% Determine the set of features that are selected by the region mask

inMask = logical(Mask(sub2ind(size(Mask), round(data.Y), round(data.X))));

% Create a new structure with only these features

data_out = structfun(@(x) x(inMask), data, 'UniformOutput', false);

% Find the distance between each pair of features

[iLoc, jLoc, dist] = position_analysis(data_out, XScale);

% Determine nearest neighbors from radial distribution function

nn = find_nearest_neighbors(iLoc, jLoc, dist, r_0, show_images);

% Add new elements to the data structure to track number of nearest
% neighbors and those neighbor's feature indices

data_out.Num_NN = sum(nn ~= 0);

for i = 1:length(nn)
    nn_pts = (nn(i, :) ~= 0); %#ok<*AGROW>
    data_out.NN_Index{i} = data_out.Index(nn_pts);
end

    