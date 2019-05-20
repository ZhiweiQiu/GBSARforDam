function c = enscoh(P)
%ENSCOH   Ensemble coherence
%   C = ENSCOH(P) returns the absolute value of the 
%   ensemble coherence given N phase residuals in
%   vector P.  The ensemble coherence is defined as
%     c = abs(sum(complex(cos(P),sin(P)))) ./ length(P);
%   If P is complex, returned is: 
%     c = abs(sum(P)) ./ length(P);
%
%   This function is used for informational purposes only,
%   it is not used during computations.
%
%   See also STUN.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% See also Fig. 2.3 in the book.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(1,1,nargin))
if (size(P,1)~=1 & size(P,2)~=1)
  error('P must be a vector');
end


% --- Compute ensemble coherence ----------------------------
if (isreal(P)) 
  c = abs(sum(complex(cos(P),sin(P)))) ./ length(P);
else
  c = abs(sum(P)) ./ length(P);
end
%%% EOF
