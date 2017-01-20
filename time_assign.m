function fitresult_out = time_assign(fitresult, XScale, Period, YPixels, SlowScan, FastScan)

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
% time_assign assigns time values to the features in the fitresult structure
% returned by find_molecules. Time is time since the beginning of the image,
% and Fast_Time is time since the beginning of the fast scan line. 
% time_assign does this by assuming that each pixel taxes a constant amount 
% of time to acquire (in seconds), reported by the scan controller as 
% Period. If a time basis was measured for the STM image when it was taken, 
% feature times should be determined from that instead, which will be more
% accurate.
%
%     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
%     DHCT.m
%         time_assign.m (this script)
%
%     Dependencies
%         none
% 
%     Dependents (*called by) 
%         DHCT.m (*)
%

% Condition input data and prepare useful constants 

fitresult_out = fitresult;

if ~iscell(fitresult_out)
    fitresult_out = {fitresult_out};
end

nm_to_m = 1e-9;

scanspeed = XScale / Period / nm_to_m;
fast_linetime = Period * YPixels;
slowscan = scanspeed / YPixels / 2;
image_time = Period * 2 * YPixels.^2;

% For loop handles trace and retrace image data, if applicable

for i = 1:length(FastScan)
    
    % Round feature locations to the nearest pixel and convert to nanometers.
    
    if SlowScan == 2 || SlowScan == 3
        fast_point = round(fitresult_out{i}.X) * XScale / nm_to_m;
        slow_point = round(fitresult_out{i}.Y) * XScale / nm_to_m;
    else
        fast_point = round(fitresult_out{i}.Y) * XScale / nm_to_m;
        slow_point = round(fitresult_out{i}.X) * XScale / nm_to_m;
    end
    
    % Determine fast-scan time from fast-scan coordinate
    
    if FastScan(i) == 1 || FastScan(i) == 2
        fast_time = fast_linetime - (fast_point / scanspeed);
    else
        fast_time = fast_point / scanspeed;
    end

    % Determine slow-scan time from slow-scan coordinate
    
    if SlowScan == 1 || SlowScan == 2
        slow_time = image_time - (slow_point / slowscan);
    else
        slow_time = slow_point / slowscan;
    end
    
    % The retrace image is delayed from the trace image by one fast-scan
    % line acquisition time
    
    if i == 2
        delta = fast_linetime;
    else
        delta = 0;
    end
    
    % Build total time from the three times determined above
    
    t_point = fast_time + slow_time + delta;
    
    % Assign time to output structure
    
    fitresult_out{i}.Time = t_point;
    fitresult_out{i}.Fast_Time = fast_time;
end

% If data was conditioned, return it to its original style.

if ~iscell(fitresult)
    fitresult_out = fitresult_out{1};
end
