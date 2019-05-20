function [X,Y] = simpos(NPS, SIZE)
%SIMPOS   Simulate position of PS points.
%   [X,Y] = SIMPOS(N) returns lying vectors with unique
%   positions in meters in area of 10 x 10 km^2.
%
%   [X,Y] = SIMPOS(N,SIZE) uses an area of SIZE(1) x SIZE(2) m^2.
%   If SIZE is a scalar a square area is simulated.
%
%   Example:
%     [X,Y] = simpos(200);
%     plot(X,Y,'r+');
%     title('PS spatial distribution');
%     xlabel('ground-range [m]');
%     ylabel('azimuth [m]');
%
%   See also STUN, SIMACQ, SIMPHI, ILSDEMO2D.

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


% --- Check input -----------------------------------
error(nargchk(1,2,nargin))
if (nargin==1)
  SIZE = [10000, 10000];% default nY,nX
end
if (length(SIZE)~=2)
  SIZE = [SIZE(1), SIZE(1)];
end
SIZE = SIZE-1;% output in [1:SIZE]



% -- Make unique, sorted positions ------------------
YX = [];
rand('state',sum(100*clock));% change state of random number generator
while (size(YX,1)~=NPS)
  N  = NPS-size(YX,1);
  YX = [YX; 1+round([SIZE(2)*rand(N,1), SIZE(1)*rand(N,1)])];% in [1:SIZE]
  YX = unique(YX,'rows');
end
Y = YX(:,1).';% return lying vector
X = YX(:,2).';% return lying vector

%%%EOF
