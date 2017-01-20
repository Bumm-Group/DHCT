function [soln, driftX, driftY, scaleFactor] = molecule_local_distortion(...
    fitresult, XScale, Period, YPixels, spacing, feat_pos_nn, varargin)

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
%     Homer L. Dodge Depaetment of Physics & Astronomy
%     The University of Oklahoma
%     Norman, OK 73019
%     bumm@ou.edu
%
% SUMMARY
% molecule_local_distortion measures local distortion at each surface
% feature's location using its nearest neighbors by approximating their
% position as a linear transformation of the lattice's structure
% and fitting the distortion parameters to a model for how the STM image 
% is distorted by the inaccuracies of the STM controller. These distortion 
% parameters can be used to recalculate the (x, y) coordinates of the 
% image and image features. 
%
% If method is 'Both', molecule_local_distortion determines hysteresis, creep
% and thermal drift simultaneously. If method is 'Hyst',
% molecule_local_distortion determines only hysteresis. If method is 'Drift',
% molecule_local_distortion determines creep and thermal drift.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         molecule_local_dstortion.m (this script)
%             drift_models.m
%
%     Dependencies
%         drift_models.m
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

% Orgnaize input variables

if nargin > 6
    if nargin > 7
        method = varargin{1};
        show_images = varargin{2};
    elseif ~ischar(varargin{1})
        method = varargin{1};
        show_images = 0;
    else
        method = 'Both';
        show_images = varargin{1};
    end
else
    method = 'Both';
    show_images = 0;
end

if iscell(fitresult)
    NumImages = length(fitresult);
else
    fitresult = {fitresult};
    NumImages = 1;
end

% Define some useful constants

nm_to_m = 1e-9;

r_0 = spacing / nm_to_m;

scanspeed = XScale / Period / nm_to_m;
slowscan = scanspeed / YPixels / 2;

options = optimoptions('lsqnonlin', 'Display', 'off');

for i = 1:NumImages  

    cen_x = feat_pos_nn{i}.cen_x;
    cen_y = feat_pos_nn{i}.cen_y;
    nn_x = feat_pos_nn{i}.nn_x;
    nn_y = feat_pos_nn{i}.nn_y;
    feat_set = feat_pos_nn{i}.feat_set;
    
    % If the drift results are problematic, we will mark that point as an
    % outlier
    
    outlier = zeros(length(feat_set), 1);
    
    % Setup variables for the parallel loop
    
    tPts{i}(1, :) = fitresult{i}.Time(feat_set); %#ok<*AGROW>
    tPts_h{i}(1, :) = fitresult{i}.Fast_Time(feat_set);
    
    driftX_t = [];
    driftY_t = [];
    scaleFactor_t = [];
    
    L = length(feat_set);
    
    % Initialize parallel pool
    
    gcp;
    
    fprintf('Progress:\n');
    fprintf([repmat(' ', 1, 20) '25%% v' repmat(' ', 1, 20) '50%% v' ...
        repmat(' ', 1, 20) '75%% v' repmat(' ', 1, 19) '100%% v' '\n\n']);
     
    parfor j = 1:L
        
        % Ellipse equation and guess for fitting to NNs
        
        fun = @(s) ((nn_x(:, j) - s(4)).*cos(s(3)) - (nn_y(:, j) - s(5))*sin(s(3))).^2/s(1).^2 ...
            + ((nn_x(:, j) - s(4)).*sin(s(3)) + (nn_y(:, j) - s(5))*cos(s(3))).^2/s(2).^2 - 1;
        guess = [r_0 r_0 0 cen_x(j) cen_y(j)];
        lb = [0, 0, -pi/2 + eps, cen_x(j) - r_0, cen_y(j) - r_0];
        ub = [Inf, Inf, pi/2, cen_x(j) + r_0, cen_y(j) + r_0];
        
        % Determine best-fit ellipse to NNs
        
        soln = lsqnonlin(fun, guess, lb, ub, options);
        
        % Give best-fit ellipse parameters reasonable variable names
        
        majorAxis = max(soln(1:2));
        minorAxis = min(soln(1:2));
        if majorAxis == soln(1)
            theta = mod(soln(3) + pi/2, pi) - pi/2;
        else % if majorAxis == soln(2)
            theta = mod(soln(3), pi) - pi/2;
        end
        
        % Change best-fit ellipse parameters into a more useful set of
        % parameters for determining the transformation
        
        height = sqrt((majorAxis * sin(theta)).^2 + (minorAxis * cos(theta)).^2);
        width = majorAxis * minorAxis / height;
        shear = height * sin(theta) * cos(theta) * (majorAxis.^2 - minorAxis.^2) / ...
            ((majorAxis * sin(theta)).^2 + (minorAxis * cos(theta)).^2);
        
        % Using those parameters, get transformation matrix elements
        
        scaleFactor_t(j) = width / r_0;
        driftX_t(j) = shear / height;
        driftY_t(j) = (r_0 - height) / height;
        
        % If the drift measurement was higher than the slow scan speed,
        % assume measurement was an outlier
        
        if or(abs(driftX_t(j)) > 1, abs(driftY_t(j)) > 1)
            outlier(j) = 1;
        end
        
        if floor(mod(j, L / 100)) == 0
            fprintf('\b%%\n');
        end
    end 
    
    % Discard outliers, and change drift matrix elements into velocity

    driftX{i} = -1 * driftX_t(outlier ~= 1) * slowscan;
    driftY{i} = -1 * driftY_t(outlier ~= 1) * slowscan;
    tPts{i} = tPts{i}(outlier ~= 1);
    tPts_h{i} = tPts_h{i}(outlier ~= 1);
    scaleFactor{i} = scaleFactor_t(outlier ~= 1);
end

if  ~strcmp(method, 'Drift') % 'Hyst' or 'Both'

    % Collect scale factor and time data from all images to fit together
    
    sf = [];
    th = [];
    
    for i = 1:length(scaleFactor)
        sf = [sf, scaleFactor{i}];
        th = [th, tPts_h{i}];
    end
    
    % Equation of line of best fit to hysteresis
    
    [guess, lb, ub] = drift_models('Hyst');
    [fun_hyst, fun_plot_hyst] = drift_models('FunHyst', th, sf);
    
    % Find best-fit curve to hysteresis
    
    Hyst_soln = lsqnonlin(fun_hyst, guess, lb, ub, options);
    
    % Plot trace and retrace scale factor measurements and best-fit curve 
    
    if show_images == 1
        Color = {'Red', 'Blue'};
        
        figure;
        hold on;
        for i = 1:length(scaleFactor)
            scatter(tPts_h{i}, scaleFactor{i}, 1, Color{i}, '.');
        end
        
        fit_hyst = fun_plot_hyst(Hyst_soln);
        [~, order] = sort(th);
        
        plot(th(order), fit_hyst(order), 'Color', 'black', 'LineWidth', 3);
        drawnow;
    end
end

if ~strcmp(method, 'Hyst') % 'Drift' or 'Both'
    
    % Collect drift and time data from all images to fit together
    
    dx = [];
    dy = [];
    td = [];
    
    for i = 1:NumImages
        dx = [dx, driftX{i}];
        dy = [dy, driftY{i}];
        td = [td, tPts{i}];
    end
    
    % Get fitting parameters from model file.
    
    [guess, lb, ub] = drift_models('Drift', slowscan);
    [fun_drift, fun_plot_x, fun_plot_y] = drift_models('FunDrift', td, dx, dy);
    
    % Find best-fit curve to drift.
    
    Drift_soln = lsqnonlin(fun_drift, guess, lb, ub, options);
    
    % Plot trace and retrace drift measurements and best-fit curve
    
    if show_images == 1
        Color = {'Red', 'Blue'};
        
        [~, order] = sort(td);
        
        figure;
        hold on;
        for i = 1:NumImages
            scatter(tPts{i}, driftX{i}, 1, Color{i}, '.');
        end
        
        fit_x = fun_plot_x(Drift_soln);
        
        plot(td(order), fit_x(order), 'Color', 'black', 'LineWidth', 3);
        drawnow;
        
        figure;
        hold on;
        for i = 1:NumImages
            scatter(tPts{i}, driftY{i}, 1, Color{i}, '.');
        end
        
        fit_y = fun_plot_y(Drift_soln);
        
        plot(td(order), fit_y(order), 'Color', 'black', 'LineWidth', 3);
        drawnow;
    end
end

% If input data was conditioned, return it to the original state

if NumImages == 1
    driftX = driftX{1};
    driftY = driftY{1};
    scaleFactor = scaleFactor{1};
end

% Output soln depends on method 

if strcmp(method, 'Drift')
    soln = Drift_soln;
elseif strcmp(method, 'Hyst')
    soln = Hyst_soln;
else
    soln = [Drift_soln Hyst_soln];
end

