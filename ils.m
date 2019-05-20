function [zfixed,sqnorm,ierr] = ils(zfloat,Linv,Dinv,Chi2,ncands,maxloop)
%ILS   Integer ambiguity resolution search
%   ZFIXED = ILS(ZFLOAT,LINV,DINV) returns the vector with integer 
%   elements closest to the given float vector ZFLOAT in a 
%   least-squares sense.  ZFLOAT is a standing vector with 
%   the float solution.  It is best to perform this search on
%   decorrelated ambiguities, i.e., ZFLOAT=Z.'*a, where 
%   Z is the Z-transformation matrix and a the float solution
%   for the ambiguities.
%   LINV and DINV are the inverses of the LtDL-decomposition 
%   of the decorrelated vc-matrix of the ambiguities, 
%   LINV=inv(L), DINV=1./D.
%
%   ZFIXED = ILS(ZFLOAT,LINV,DINV,CHI2) uses CHI2 as bound
%   on the search space (ellipsoid).  If CHI2 is not given
%   an extended bootstrap estimation is performed to obtain
%   a bound such that at least 1 solution is contained.
%
%   ZFIXED = ILS(ZFLOAT,LINV,DINV,CHI2,NCANDS) returns NCANDS
%   solutions in a matrix, in order of ascending squared 
%   norm of the solutions.
%
%   ZFIXED = ILS(ZFLOAT,LINV,DINV,CHI2,NCANDS,MAXLOOP) stops
%   the search if MAXLOOP ambiguities have been sought.  For
%   PS processing, the default, MAXLOOP=N.^3, may be reasonable 
%   to prevent extremely long searches with noisy data. Use
%   MAXLOOP=0 to prevent premature exit.
%
%   [ZFIXED, SQNORM] = ILS(...) returns the corresponding squared
%   norm.  If NCANDS is larger than 1 SQNORM is a sorted vector.
%
%   [ZFIXED, SQNORM, IERR] = ILS(...) returns an error code
%   IERR=0 for successful exit, IERR=IERR+1 if not enough candidates 
%   could be found, IERR=IERR+2 if MAXLOOP was reached.
%
%   Example:
%      Q      = [4.25, -0.5, -1.0,  0.5; ...
%               -0.50, 21.0, -2.0, 15.0; ...
%               -1.00, -2.0,  2.0, -4.0; ...
%                0.50, 15.0, -4.0, 21.0];
%      afloat = [0.1; 0.6; 0.2; 0.4];
%
%      Compute prelimanry fixed solution using the bootstrap
%      in order to obtain a bound for the search space:
%        [Z,L,D,Qd,zfloat]         = zt(Q,afloat);
%        ncands = 2;
%        [zfixed_bs, zsqnorm_bs]   = ebs(zfloat,L,D,ncands);
%      Use this bound to search the solutions space more thoroughly:
%        Linv = inv(L);
%        Dinv = 1./D;
%        [zfixed_ils, zsqnorm_ils] = ils(zfloat,Linv,Dinv,max(zsqnorm_bs),ncands);
%
%      disp('squared norm using decorrelated ambiguities using bootstrap: ');
%        zsqnorm_bs  = zsqnorm_bs
%      disp('squared norm using decorrelated ambiguities using ILS: ');
%        zsqnorm_ils = zsqnorm_ils
%      disp('solution for decorrelated ambiguities using bootstrap:');
%        zfixed_bs   = zfixed_bs
%      disp('solution for decorrelated ambiguities using ILS:');
%        zfixed_ils  = zfixed_ils
%
%   See also STUN, EBS, LTDL, ZT, BITAND.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This function pertains to Chapter 3.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.5 $  $Date: 2005/12/07 08:23:06 $

% This routines is based on the ILS search algorithm
% in the lambda toolbox developed by Delft University of Technology,
% originally distributed as:
% ----------------------------------------------------------------------
% File.....: lsearch.m
% Date.....: 19-MAY-1999
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------
%
% Changes with respect to that version are related to using 
% matrix notation (vectorization), interface changes to improved speed
% for PS processing due to passing the inverse of L and D and 
% due to not sorting the candidates.  A maxloop counter
% is introduced to stop the search when more than the specified
% number of ambiguities are sought.  Bugs in the original code
% were related to not initializing variables, or not with
% the correct value (e.g., Inf).



% --- Check input (remove this for faster code) ----------------
error(nargchk(3,6,nargin))
% --- Check dimensions of input --------------------------------
if (diff(size(Linv))~=0)
  error('Linv must be square');
end
if (prod(size(Dinv))~=size(Linv,1))
  error('Dinv must be vector corresponding to Linv');
end
if (size(zfloat,1)~=size(Linv,1) | size(zfloat,2)~=1)
  error('zfloat must be standing vector corresponding to Linv');
end
% --- Set defaults ---------------------------------------------
if (nargin < 6)
  maxloop = length(zfloat)^3.0;% prevent extremly long searches
end
if (maxloop==0)
  maxloop = 2^31;
end
if (nargin < 5)
  ncands = 1;% default 1 candidate sought
end
if (nargin < 4)% use bootstrap estimator to obtain bound
  [zfixed_bs, zsqnorm_bs] = ebs(zfloat,inv(Linv),1./Dinv,ncands);
  Chi2 = max(zsqnorm_bs);
end



% --- Prevent not finding a solution at boundary due to --------
% --- numerical rounding, e.g., if Chi2 was computed using ebs -
Chi2 = Chi2 + 1d-6;% do not use Matlab's "eps" (too small)
% --- Initialization -------------------------------------------
n         = length(zfloat);
right     = [zeros(n,1); Chi2];
left      = [zeros(n+1,1)];
dq        = [Dinv(2:n)./Dinv(1:n-1), 1./Dinv(n)];



% --- Loop counters --------------------------------------------
ncan      = 0;
i         = n + 1;
iold      = i;
ierr      = 0;% status
cntloop   = 0;% exit search when it takes too long

% --- Guarantee boolean support also for Matlab 5.3 ------------
% --- In later version "true" is a built-in function -----------
% --- Note that I also use '&' not '&&' in order to be ---------
% --- compatible with Matlab 5.3 -------------------------------
true  = logical(1); 
false = logical(0);
% --- Loop booleans --------------------------------------------
cand_n    = false;
c_stop    = false;
endsearch = false;



% --- Allocate and initialize matrices -------------------------
zfixed    = zeros(n,ncands);
% --- It is important to initialize sqnorm with a large 
% --- number in case not enough candidates are found.
sqnorm    = Inf+zeros(1,ncands);% set to infinite
% --- Initializing dist prevents a bug in original code, i.e.,
% --- usage of nn toolbox::dist() function if it exists.
dist      = zeros(n,1);
endd      = zeros(n,1);
lef       = zeros(n,1);



% --- Start search ----------------------------------------------
while (endsearch==false);
  % --- break loop if it takes too long -------------------------
  cntloop=cntloop+1;
  if (cntloop > maxloop)
    endsearch = true;
    %warning(['maximum number of loops reached > ', num2str(maxloop)]);
    ierr = ierr+2;% set second bit of error flag
  end
  %
  % --- Start the main search loop ------------------------------
  i = i-1;
  if (iold <= i)
    lef(i) = lef(i) + Linv(i+1,i);
  else 
    if (i+1 <= n)
      lef(i) = Linv(i+1:n,i).'*dist(i+1:n);
    end
  end
  iold      = i;
  right(i)  = (right(i+1) - left(i+1)) * dq(i);
  reach     = sqrt(right(i));
  delta     = zfloat(i) - reach - lef(i);
  dist(i,1) = ceil(delta) - zfloat(i);
  if (dist(i,1) > reach-lef(i))
    % --- There is nothing at this level, so backtrack ----------
    cand_n = false;
    c_stop = false;
    while ((c_stop==false) & (i<n))
      i = i + 1;
      if (dist(i) < endd(i))
        dist(i) = dist(i) + 1;
        left(i) = (dist(i) + lef(i)) ^ 2;
        c_stop = true;
        if (i == n)
	  cand_n = true;
	end
      end
    end
    if ((i==n) & (cand_n==false))
      endsearch = true;
    end
  else
    % --- Set the right border ---------------------------------
    endd(i) = reach - lef(i) - 1;
    left(i) = (dist(i,1) + lef(i))^2;
  end

  if (i == 1)
    % --- Collect the integer vectors and corresponding norms --
    % --- add to vectors "zfixed" and "sqnorm" if --------------
    % --- * Less then "ncands" candidates found so far ---------
    % --- * The squared norm .lt. a_previous_one ---------------
    t       = Chi2 - (right(1)-left(1)) * Dinv(1);
    endd(1) = endd(1) + 1;
    while (dist(1) <= endd(1))
      if (ncan < ncands)
        ncan             = ncan + 1;
        zfixed(1:n,ncan) = dist + zfloat;
        sqnorm(ncan)     = t;
      else
        [maxnorm,ipos] = max(sqnorm);
        if (t < maxnorm)
          zfixed(1:n,ipos) = dist + zfloat;
          sqnorm(ipos)     = t;
        end
      end
      t       = t + (2 * (dist(1) + lef(1)) + 1) * Dinv(1);
      dist(1) = dist(1) + 1;
    end
    % --- And backtrack ---------------------------------------
    cand_n = false;
    c_stop = false;
    while ((c_stop==false) & (i<n))
      i = i + 1;
      if (dist(i) < endd(i))
        dist(i) =  dist(i) + 1;
        left(i) = (dist(i) + lef(i))^2;
        c_stop = true;
        if (i == n)
	  cand_n = true;
	end
      end
    end
    if ((i==n) & (cand_n==false))
      endsearch = true;
    end
  end
end



% --- Sort the resulting candidates, according to the norm ---
zfixed = round(zfixed);% remove small numerical inaccuracies
if (ncands > 1)
  [sqnorm,idx] = sort(sqnorm);
  zfixed       = zfixed(:,idx);
end



% --- Check for errors ---------------------------------------
if (ncan < ncands)
  %warning(['Not enough candidates.  Only found ', num2str(ncan)]);
  ierr = ierr+1;% set first bit of error flag
end

%%%EOF
