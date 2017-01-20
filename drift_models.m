function varargout = drift_models(mode, varargin)

%     ver 12/11/2016
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
% drift_models contains the fitting functions to hysteresis, creep, and
% thermal drift used by a variety of functions in DHCT. They have been
% collected in one place for ease of modification and checking for self-
% consistency. 'Symfun' mode returns a symbolic function for integration
% that is used in fitresult_modify and ImCoords. The other four modes are
% used for fitting best-fit curves to the data in molecule_local_distortion.
% "FunHyst' and 'FunDrift' contain inline functions for fitting and
% plotting the best-fit curves, while 'Hyst' and 'Drift' contain the guess,
% lower bound, and upper bound for fitting with lsqnonlin.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         fitresult_modify.m
%            -drift_models.m (this script)
%         ImCoords.m
%            -drift_models.m (this script)
%         molecule_local_distortion.m
%            -drift_models.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         fitresult_modify.m (*)
%         ImCoords.m (*)
%         molecule_local_distortion.m (*)
% 

if strcmpi(mode, 'Hyst')
% [guess, lb, ub] = drift_models('Hyst');
    
    % Guess, and lower and upper bounds, for fitting function.
    
    guess = [1 0.1];
    lb = [-Inf -Inf];
    ub = [Inf Inf];
    
    varargout{1} = guess;
    varargout{2} = lb;
    varargout{3} = ub;
    
elseif strcmpi(mode, 'FunHyst') 
% [fun_hyst, fun_plot_hyst] = drift_models('FunHyst', time_h, scaleFactor);
    if nargin == 3
        
        time_h = varargin{1};
        scaleFactor = varargin{2};
        
        % Fitting functions for hysteresis model
        
        fun_hyst = @(s) s(1) * time_h .^ s(2) - scaleFactor;
        fun_plot_hyst = @(s) s(1) * time_h .^ s(2);
        
        varargout{1} = fun_hyst;
        varargout{2} = fun_plot_hyst;
    else
        error('Invalid number of input arguments for "FunHyst" mode.');
    end
    
elseif strcmpi(mode, 'Drift')
% [guess, lb, ub] = drift_models('Drift', slowscan);
    if nargin == 2
        
        slowscan = varargin{1};
        
        % Guess, and lower and upper bounds, for fitting function.
        
        guess = [0 0 0 0 0 0 10];
        lb = [-slowscan -slowscan -10 -slowscan -slowscan -10 0];
        ub = [slowscan slowscan 10 slowscan slowscan 10 1000];
        
        varargout{1} = guess;
        varargout{2} = lb;
        varargout{3} = ub;
        
    else
        error('Invalid number of input arguments for "Drift" mode.');
    end
    
elseif strcmpi(mode, 'FunDrift')
% [fun_drift, fun_plot_x, fun_plot_y] = drift_models('FunDrift', td, dx, dy);
    if nargin == 4
        
        time_d = varargin{1};
        driftX = varargin{2};
        driftY = varargin{3};
        
        % Fitting functions for drift model
        
        fun_drift = @(s) [s(1) + s(2) * time_d + s(3) ./ (s(7) + time_d) - driftX, ...
            s(4) + s(5) * time_d + s(6) ./ (s(7) + time_d) - driftY];
        fun_plot_x = @(s) s(1) + s(2) * time_d + s(3) ./ (s(7) + time_d);
        fun_plot_y = @(s) s(4) + s(5) * time_d + s(6) ./ (s(7) + time_d);
        
        varargout{1} = fun_drift;
        varargout{2} = fun_plot_x;
        varargout{3} = fun_plot_y;
    else
        error('Invalid number of input arguments for "Fundrift" mode.');
    end
    
elseif strcmp(mode, 'Symfun')
% [xDrift(t), yDrift(t), h(t)] = drift_models('Symfun', solution_vector);
    if nargin == 2
               
        % Coefficients of drift model in X and Y
        
        Xlindrift = varargin{1}(1);
        Xacc = varargin{1}(2);
        Xamp = varargin{1}(3);
        Ylindrift = varargin{1}(4);
        Yacc = varargin{1}(5);
        Yamp = varargin{1}(6);
        t_0 = varargin{1}(7);
        
        % Coefficients of hysteresis model, numbering starts after
        % the last drift coefficient
        
        amplitude = varargin{1}(8);
        exponent = varargin{1}(9);
        
        % Define symbolic functions to integrate that are consistent with
        % the inline functions defined above
        
        syms t
        
        xDrift(t) = Xlindrift + Xacc * t + Xamp / (t + t_0);
        yDrift(t) = Ylindrift + Yacc * t + Yamp / (t + t_0);
        hyst(t) = amplitude * t .^ exponent - 1;
        
        varargout{1} = xDrift;
        varargout{2} = yDrift;
        varargout{3} = hyst;
    else
        error('Invalid number of input arguments for "Symfun" mode.');
    end
else
    error('Invalid mode. Check drift_models.m for valid modes.');
end
