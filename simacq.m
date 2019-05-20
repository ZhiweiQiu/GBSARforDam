function [bperp, btemp, acqtimes] = simacq(NIFG, SD, NY)
%SIMACQ   Simulate acquisition baselines.
%   BPERP = SIMACQ(N) returns standing vector with 
%   perpendicular baseline [meter] for N interferograms with 
%   a single master.  Baseline simulated with a standard
%   normal distribution with 400 meter standard deviation.
%
%   BPERP = SIMACQ(N, SD) uses SD as standard deviation to
%   generate the baselines.
%
%   [BPERP, BTEMP] = SIMACQ(N) additionally returns the
%   temporal baseline [year].  Orbit revolution of 35 days
%   is assumed.  Data simulated over 8 years time span
%   with an uniform distributon.  BTEMP is sorted.
%
%   [BPERP, BTEMP] = SIMACQ(N, SD, NY) uses NY time span.
%
%   [BPERP, BTEMP, ACQ] = SIMACQ(...) additionally returns
%   the acquisition times, as fractional year, e.g., 1996.4,
%   in a standing vector of length N+1, where the master
%   image is the first element.  The master is the median
%   acquisition.
%
%   Example:
%     [bperp, btemp, acqt] = simacq(30);
%     plot([0;bperp],acqt,'r+');% add master
%     title('Baseline distribution');
%     xlabel('perpendicular baseline');
%     ylabel('acquisition time');
%
%   See also STUN, SIMPHI, SIMPOS, ILSDEMO2D.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $


% --- Check input --------------------------------------
error(nargchk(1,3,nargin))
if (nargin<=1)
  SD = 400;% default
end
if (nargin<=2)
  NY = 8;% default
end


% --- Simulate the baselines ---------------------------
bperp = simbperp(NIFG, SD);
[btemp, acqtimes] = simbtemp(NIFG, NY);




% --- Subfunction simbperp -----------------------------
% --- return perpendicular baselines as normally -------
% --- distributed --------------------------------------
function bperp = simbperp(NIFG,SD)
error(nargchk(2,2,nargin))
randn('state',sum(100*clock));% change state of random number generator
bperp = SD*randn(NIFG,1);% w.r.t. master
%EOF



% --- Subfunction simbtemp -----------------------------
% --- return temporal baselines and acquisition times --
% --- as fractional year.  Times are uniform over ------
% --- NYEAR using 35 day repeat orbit;  The master is --
% --- the median acq. time, returned as first element --
% --- in vector ACQ_TIMES of length NIFG+1 -------------
% --- Output further time sorted ascending -------------
% --- It is not guaranteed that btemp is unique --------
function [btemp, acq_times] = simbtemp(NIFG,NYEAR)
error(nargchk(2,2,nargin))
orbit_cycle      = 35/365;% [year] repeat interval
first_acq        = 1992.5;% first acquisition
rand('state',sum(100*clock));% change state of random number generator
acq_times        = sort(rand(NIFG+1,1));% in [0,1]
acq_times        = NYEAR*(acq_times-min(acq_times))/ ...
                     (max(acq_times)-min(acq_times));% [0,NYEAR]
acq_times        = first_acq + round(acq_times/orbit_cycle)*orbit_cycle;%
idx              = find(diff(acq_times)==0);% find same dates
acq_times(idx)   = acq_times(idx)-1/365;% simulate tandem image
idx_m            = round(NIFG/2);
acq_m            = acq_times(idx_m);%master central in time
acq_times(idx_m) = [];% remove master temporarily
acq_times        = [acq_m; acq_times];% add master as first
btemp            = acq_times-acq_m;
btemp(1)         = [];% remove master from temporal baseline
%EOF



%%%EOF
