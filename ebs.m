function [zfixed, sqnorm] = ebs(zfloat,L,D,ncands, INVLTDL)
%EBS   Extended bootstrap estimator for ambiguity resolution.
%   ZFIXED = EBS(ZFLOAT,L,D) returns the bootstrap
%   estimator for the fixed ambiguities based on the float 
%   solution in standing vector ZFLOAT and the decomposition 
%   L, D.  For best results, zfloat. L, and D should be 
%   Z-transformed, i.e., decorrelated.
%
%   ZFIXED = EBS(ZFLOAT,L,D,NCANDS) returns the best 
%   NCANDS fixed solutions in a matrix with NCANDS columns.
%
%   ZFIXED = EBS(ZFLOAT,L,D,NCANDS,INVQ) with INVQ=inv(L.'*diag(D)*L)
%   prevents computation of this matrix in this routine, which 
%   makes it faster if this routine is called many times with the same
%   model (Qz_hat matrix).
%
%   [ZFIXED, SQNORM] = EBS(...) optional returns the 
%   squared norm in a vector of length NCANDS.  The (largest)
%   squared norm of the bootstrap solution can be used to bound
%   the search space during the ILS search.
%
%   The extended bootstrap method uses a series of slightly
%   altered bootstraps to compute a series of fixed solutions.
%   For N float ambiguities, N+2 solutions are computed, and
%   the best NCANDS solutions are returned with the smallest norm.
%
%   The bootstrap solution is faster than the integer least-squares
%   search, but is not guaranteed to find the minimum norm solution.
%
%   Example:
%      Q      = [4.25, -0.5, -1.0,  0.5; ...
%               -0.50, 21.0, -2.0, 15.0; ...
%               -1.00, -2.0,  2.0, -4.0; ...
%                0.50, 15.0, -4.0, 21.0];
%      afloat = [0.1; 0.6; 0.2; 0.4];
%
%      Best fixed solution using the original ambiguities:
%        [L,D] = ltdl(Q);
%        [afixed, asqnorm] = ebs(afloat,L,D,1);
%      disp('squared norm using original ambiguities should be 0.0319');
%        asqnorm = asqnorm
%      disp('solution for fixed ambiguities should be [0 1 0 1]:');
%        afixed  = afixed
%
%      Best fixed solution using the decorrelated ambiguities:
%        [Z,L,D,Qd,zfloat] = zt(Q,afloat);% L,D are Z-transformed too
%        [zfixed, zsqnorm] = ebs(zfloat,L,D,1);
%      disp('squared norm using decorrelated ambiguities:');
%        zsqnorm = zsqnorm;
%      disp('solution for decorrelated ambiguities:');
%        zfixed = zfixed
%      disp('solution for fixed ambiguities should be [0 1 0 1]:');
%        afixed_Z = inv(Z.')*zfixed
%
%   See also STUN, ILS, LTDL, ZT.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% See section 3.3.1. in this book.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $

% This routines is based on the simple bootstrap algorithm
% in the lambda toolbox developed by Delft University of Technology,
% distributed as:
% File.....: chistart.m
% Date.....: 19-MAY-1999
% Modified.: 05-MAR-2001, by P. Joosten
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
%
% Changes with respect to that version are related to the 
% extention using the series of slightly altered bootstrap 
% solutions, the help, interface, speed, and a small bugfix.



% --- Check input (remove this for faster code) ----------------
error(nargchk(3,5,nargin));
n = length(zfloat);
if (nargin < 5)
  % --- Code can be sped-up by offering INVLTDL to this routine -----
  % --- which is inv(Q);  This save time  ---------------------------
  % --- in case this routine is called many times with same input ---
  INVLTDL    = inv(L.'*diag(D)*L);% or precompute outside loop
end
if (nargin < 4)
  ncands = 1;% default
end
if (diff(size(L)) ~= 0)
  error('L must be square');
end
if (prod(size(D)) ~= size(L,1))
  error('D must be vector corresponding to L');
end
if ((size(zfloat,1)~=size(L,1)) | (size(zfloat,2) ~= 1))
  error('zfloat must be a standing vector');
end
if (ncands > n)
  error('too many candidates requested');
end





% --- Initialize tmp variables to store solution --------------
all_chi    = zeros(1,n+2);% squared norms of all solutions
all_zfixed = zeros(n,n+2);% all solutions for each altered bootstrap



% --- Perform the series of altered bootstraps ----------------
% --- k loops over n+2 slightly altered series of bootstraps.
for k = n+2:-1:1
  this_zfloat = zfloat;% reset to original float solution
  this_zfixed = zeros(n,1);% standing vector, solutions for this bootstrap
  % --- D is smallest for first element -----------------------
  % --- so it would seem better to start with 1 ---------------
  start_idx    = n;% However, start with the last one and alter it
  bootstrapped = round(this_zfloat(start_idx));% 
  % --- Adapt the first bootstrapped estimate first two times ----------
  switch k
    case n+2,   this_zfixed(start_idx)=bootstrapped+1;
    case n+1,   this_zfixed(start_idx)=bootstrapped-1;
    otherwise,  this_zfixed(start_idx)=bootstrapped;
  end
  % --- Now bootstrap the others based on first ------------------------
  % --- i is counter of next bootstrapped ambiguity --------------------
  for i = n-1:-1:1
    % --- a is standing vector; elements (i+1:n) are already fixed ---
    dw             = L(i+1:n,i).' * (this_zfloat(i+1:n) - this_zfixed(i+1:n));
    this_zfloat(i) = this_zfloat(i) - dw;
    bootstrapped   = round(this_zfloat(i));% estimate for ambiguity i
    % --- Change this bootstrapped solution to little worse ------------
    if (i ~= k)
      this_zfixed(i) = bootstrapped;% normal case
    else
      if (this_zfloat(i) > zfloat(i))
        this_zfixed(i) = bootstrapped + 1;% adapt this series
      else
        this_zfixed(i) = bootstrapped - 1;% adapt this series
      end
    end
  end
  % --- Store this solution and the squared norm ----------------
  tmp             = zfloat-this_zfixed;
  all_chi(k)      = tmp.' * INVLTDL * tmp;
  all_zfixed(:,k) = this_zfixed;
end



% --- Return the best ncands solutions and norms ---------------------
if (ncands == 1)
  [sqnorm, idx]  = min(all_chi);
  zfixed         = all_zfixed(:,idx);
else
  [all_chi, idx] = sort(all_chi);
  sqnorm         = all_chi(1:ncands);
  zfixed         = all_zfixed(:,idx(1:ncands));
end
 
%%%EOF

