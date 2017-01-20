function fitresult = fmgaussfit_reduced(xx, yy, zz, initial_values)
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface.
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% link to the mathworks download page http://www.mathworks.com/matlabcentral/fileexchange/41938-fit-2d-gaussian-with-optimization-toolbox

% Edited 2016, Mitchell Yothers
%
% Removed parts of this code that were unneccessary for molecule finding to
% increase speed and modified assumptions to more accurately reflect 
% assumptions for STM images. 

%% Condition the data
xData = xx;
yData = yy;
xyData = {xData, yData};
zData = zz;

%% Set up the startpoint 

amp = max(zData) - min(zData); % amp is the amplitude. 

if initial_values == 0
    ang = 0; % angle in degrees.
    sx = 5;
    sy = 5;
    x0 = mean(xData); % guess that it is at the mean
    y0 = mean(yData); % guess that it is at the mean
else
    ang=initial_values(2);
    sx=initial_values(3) + 2;
    sy=initial_values(4) + 2;
    x0=initial_values(5);
    y0=initial_values(6);
end

z0 = min(zData);
xmax = max(xData) + 2; 
ymax = max(yData) + 2; 
zmax = max(zData) + amp; % amp is the amplitude.
xmin = min(xData) - 2; 
ymin = min(yData) - 2; 
zmin = min(zData) - amp;

%% Set up fittype and options. 
Lower = [0, 0, 0, 0, xmin, ymin, zmin]; 
Upper = [3 * amp, 180-eps, Inf, Inf, xmax, ymax, zmax]; % angles greater than 180 are redundant 
StartPoint = [amp, ang, sx, sy, x0, y0, z0];

tols = 1e-16;
options = optimset('Algorithm','trust-region-reflective',...
    'Display','off',...
    'MaxFunEvals',5e2,...
    'MaxIter',5e2,...
    'TolX',tols,...
    'TolFun',tols,...
    'TolCon',tols ,...
    'UseParallel','always');

fitresult = lsqcurvefit(@gaussian2D, StartPoint, xyData, zData, Lower, Upper, options);

end

function z = gaussian2D(s, xy)

% compute 2D gaussian
xdata = (xy{1} - s(5)) .* cosd(s(2)) - (xy{2} - s(6)) .* sind(s(2));
ydata = (xy{1} - s(5)) .* sind(s(2)) + (xy{2} - s(6)) .* cosd(s(2));
z = s(1) .* gaussmf(xdata, [s(3) 0]) .* gaussmf(ydata, [s(4) 0]) + s(7);

end



