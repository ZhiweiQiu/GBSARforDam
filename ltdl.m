function [L,D] = ltdl(Q)
%LTDL   LtDL-decompostion of square matrix.
%   [L,D] = LTDL(Q) returns the LtDL decomposition 
%   of the square and symmetric matrix Q, such that 
%   Q = L.'*diag(D)*L.
%
%   On output, L is the n by n factor matrix (strictly
%   lower triangular) and D is a lying vector containing
%   diagonal elements.
%
%   Example:
%      Q = [4.25, -0.5, -1.0,  0.5; ...
%          -0.50, 21.0, -2.0, 15.0; ...
%          -1.00, -2.0,  2.0, -4.0; ...
%           0.50, 15.0, -4.0, 21.0];
%      [L,D]   = ltdl(Q);
%      Q_check = L.'*diag(D)*L
%      status  = max(max(abs(Q-Q_check)));
%      if (status~=0) error('something went wrong!'); end
%
%   See also STUN, EBS, ILS, ZT.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% See section 3.1

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $

% This routines is based on the decomposition algorithm
% in the lambda toolbox developed by Delft University of Technology,
% distributed as:
% File.....: ldldecom
% Date.....: 19-MAY-1999
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
%
% Changes with respect to that version are mainly related 
% to the help and formating.


% --- Check input ------------------------------------------
error(nargchk(1,1,nargin))
if (diff(size(Q))~=0)
  error ('Input matrix must be square');
end


% --- Perform factorization --------------------------------
for i=size(Q,1):-1:1
  D(i) = Q(i,i);
  if (D(i)<0)
    error ('Input matrix is not positive definite!');
  end
  L(i,1:i) = Q(i,1:i)/sqrt(Q(i,i));
  for j=1:i-1
    Q(j,1:j) = Q(j,1:j)-L(i,1:j)*L(i,j);
  end
  L(i,1:i) = L(i,1:i)/L(i,i);
end

%%%EOF
