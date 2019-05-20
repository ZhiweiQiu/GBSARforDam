function [W,NW] = wrap(P)
%WRAP   wrap phase values to principal interval.
%   W = WRAP(P) returns the wrapped phase values in the
%   interval [-pi,pi).  If P is a matrix, so is W.
%
%   [W, NW] = WRAP(P) optionally returns the number of
%   wraps required to wrap P.
%
%   The following relations hold, where wrapped in [-pi,pi):
%
%     unwrapped = wrap(unwrapped) + numwraps(unwrapped)*2*pi;
%
%     wrapped   = unwrapped - numwraps(unwrapped)*2*pi;
%
%     wrapped = atan2(sin(P),cos(P));% alternative method
%
%   Example:
%     uw       = -10:0.1:10;% unwrapped phase signal
%     x        = 1:length(uw);
%     [w,nw]   = wrap(uw);
%     uw_check = w + nw * 2 *pi;
%     plot(x,uw,'r', x,uw_check,'r.', x,w, 'b', x,nw,'g');
%     title('r: unwrapped;  b: wrapped;  g: number of wraps');
%
%   See also STUN.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% See footnote in Section 2.2.1

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(1,1,nargin))
if (~isreal(P)) 
  error('input matrix must be real');
end



% --- Compute wrapped phase ---------------------------------
% --- too slow: W = atan2(sin(P),cos(P));
if (nargout==1)
  W  = mod(P+pi,2*pi) - pi;% use mod, not rem!
else
  NW = round(P/(2*pi));% number of wraps
  W  = P-NW*2*pi;
end



%%% EOF
