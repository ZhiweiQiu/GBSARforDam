function p = bs_success(D)
%BS_SUCCESS   Success rate for simple bootstrap estimator.
%   P = BS_SUCCESS(D) returns the theoretical success rate
%   for the simple bootstrap estimator.  This can be considered
%   a lower bound for the success rate for the extended 
%   bootstrap estimator and the ILS estimator.
%
%   Example:
%      Q = [4.25, -0.5, -1.0,  0.5; ...
%          -0.50, 21.0, -2.0, 15.0; ...
%          -1.00, -2.0,  2.0, -4.0; ...
%           0.50, 15.0, -4.0, 21.0];
%      [L,D]   = ltdl(Q);
%      srate   = bs_success(D);
%      disp(['success rate for simple bootstrap: ', num2str(srate)]);
%
%   See also STUN, NORMCDF

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



% --- Check input ------------------------------------------
error(nargchk(1,1,nargin))
if (diff(size(Q))~=0)
  error ('Input matrix must be square');
end


% --- Evaluate (Eq.3.20) -----------------------------------
x = 1./(2.*sqrt(D));
if (exist('normcdf','file'))
  p = prod(2*normcdf(x,0,1) - 1);
else
  p = prod(erfc(-x./sqrt(2))-1);
end

%%%EOF
