function nn = find_nearest_neighbors(iLoc, jLoc, dist, r_0, varargin)

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
% find_nearest_neighbors determines which image features are nearest
% neighbors from their distance distribution, and returns a sparse matrix
% whose nonzero elements contain the distance between a feature and one of
% its nearest neighbors, in the (iLoc, jLoc) position in the marix. iLoc
% and jLoc are two different feature indexes.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         data_mask_select.m
%            -find_nearest_neighbors.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         data_mask_select.m (*) 
%

if nargin > 4
    show_images = varargin{1};
else
    show_images = 0;
end

% Number of times the NN distance to select data from, for fitting and
% display

num_nn_fit = 2.25;
num_nn_display = 10;

% Generate radial distribution function from histogram of distances

r = (0.01:0.01:num_nn_fit) .* r_0;
h = hist(dist(dist <= r(end)), r) ./ r;

h = h / mean(h);

% Determine best-fit feature spacing by fitting Gaussians to NN, NNN
% and 2*NN peaks of radial distribution function. This assumes a
% trigonal lattice.

options = optimoptions('lsqnonlin', 'Display', 'off');

fun = @(s) s(1) * gaussmf(r, [s(4) s(5)]) + s(2) * gaussmf(r, [s(4) sqrt(3)*s(5)])...
    + s(3) * gaussmf(r, [s(4) 2*s(5)]) - h;
guess = [5, 5, 5, 0.05, r_0];
lb = [0, 0, 0, 0, 0.75 .* r_0];
ub = [10, 10, 10, Inf, 1.25 .* r_0];
soln = lsqnonlin(fun, guess, lb, ub, options);

% Determine NN features from threshold and best-fit NN distance

nn_threshold = [soln(5) * 0.65, soln(5) * 1.35];

if show_images
    
    % Regenerate radial distribution function from histogram of
    % distances, to plot it
    
    r = (0.01:0.01:num_nn_display) .* r_0;
    h = hist(dist(dist <= r(end)), r) ./ r;
    
    h = h / mean(h);
    
    % Plot RDF, best-fit curve, and NN threshold
    
    figure;
    bar(r, h, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    hold on;
    plot(r, soln(1) * gaussmf(r, [soln(4) soln(5)]) + soln(2) * gaussmf(r, [soln(4) sqrt(3)*soln(5)])  + soln(3) * gaussmf(r, [soln(4) 2*soln(5)]), 'LineWidth', 2, 'Color', 'Red');
    plot([nn_threshold(1), nn_threshold(1)], [max(h), 0], 'LineWidth', 2, 'Color', 'Black');
    plot([nn_threshold(2), nn_threshold(2)], [0, max(h)], 'LineWidth', 2, 'Color', 'Black');
    drawnow;
end


% Select NN from threshold and distance

nn_pairs = find(and(dist > nn_threshold(1), dist < nn_threshold(2)));
iLoc_nn(1:length(nn_pairs)) = iLoc(nn_pairs);
jLoc_nn(1:length(nn_pairs)) = jLoc(nn_pairs);
dist_nn(1:length(nn_pairs)) = dist(nn_pairs);

% Generate NN matrix

iLoc_nn_total = [iLoc_nn jLoc_nn];
jLoc_nn_total = [jLoc_nn iLoc_nn];
dist_nn_total = [dist_nn dist_nn];

nn = sparse(iLoc_nn_total, jLoc_nn_total, dist_nn_total);