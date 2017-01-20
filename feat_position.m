function feat_pos_nn = feat_position(Data, XScale, SlowScan)

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
% This script determines (x, y) coordinates of sets of image feature 
% nearest neighbors from their indices, in preparation for 
% molecule_local_distortion to determine local distortion parameters from 
% their locations.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         feat_position.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

if iscell(Data)
    num_images = length(Data);
else
    num_images = 1;
end

nm_to_m = 1e-9;

% num_nn = 6 for a trigonal lattice

num_nn = 6;

for i = 1:num_images
    
    % Determine features with correct # of NNs
    
    feat_set = find(cellfun(@(x) length(x) == num_nn, Data{i}.NN_Index));
    feat_pos_nn{i}.feat_set = feat_set;
    feat_pos_nn{i}.Index =  Data{i}.Index(feat_set);
    
    if SlowScan == 2 || SlowScan == 3
        
        % Get locations of features
        
        feat_pos_nn{i}.cen_x = Data{i}.X(feat_set) * XScale / nm_to_m;
        feat_pos_nn{i}.cen_y = Data{i}.Y(feat_set) * XScale / nm_to_m;
        for j = 1:length(feat_set)
            
            % Get indices of those features' NNs
            
            nn_pts = [];
            nn_ind = Data{i}.NN_Index{feat_set(j)};
            for k = 1:length(nn_ind)
                nn_pts(k) = find(Data{i}.Index == nn_ind(k));
            end
            
            % Get locations of NNs
            
            feat_pos_nn{i}.nn_x(:, j) = Data{i}.X(nn_pts) * XScale / nm_to_m; %#ok<*AGROW>
            feat_pos_nn{i}.nn_y(:, j) = Data{i}.Y(nn_pts) * XScale / nm_to_m;
        end
    else
        
        % Get locations of features
        
        feat_pos_nn{i}.cen_y = Data{i}.X(feat_set) * XScale / nm_to_m;
        feat_pos_nn{i}.cen_x = Data{i}.Y(feat_set) * XScale / nm_to_m;
        for j = 1:length(feat_set)
            
            % Get indices of those features' NNs
            
            nn_pts = [];
            nn_ind = Data{i}.NN_Index{feat_set(j)};
            for k = 1:length(nn_ind)
                nn_pts(k) = find(Data{i}.Index == nn_ind(k));
            end
            
            % Get locations of NNs
            
            feat_pos_nn{i}.nn_y(:, j) = Data{i}.X(nn_pts) * XScale / nm_to_m;
            feat_pos_nn{i}.nn_x(:, j) = Data{i}.Y(nn_pts) * XScale / nm_to_m;
        end
    end
end