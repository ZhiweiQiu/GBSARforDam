function [Z,L,D,ZtQZ,zfloat] = zt(Q,afloat)
%ZT   Z-transform (decorrelate) a vc-matrix of ambiguities.
%   Z = ZT(Q) returns the Z-transformation matrix to decorrelate
%   matrix Q.
%
%   [Z,L,D] = ZT(Q) additionally returns matrix L and diagonal D
%   of the LtDL-decomposition of the decorrelated Q matrix (not
%   of Q itself!).
%
%   [Z,L,D,QZ] = ZT(QA) additionally returns the decorrelated vc-matrix
%   QZ = Z.' * QA * Z.
%
%   [Z,L,D,QZ,ZFLOAT] = ZT(QA, AFLOAT) with additional input the
%   standing vector of ambiguities, additionally returns the 
%   Z-transformed vector of ambiguities ZFLOAT = Z.' * AFLOAT.
%   
%   Example:
%      Qa = [4.25,  -0.50,  -1.0,    0.5; ...
%           -0.50,  21.0,   -2.0,   15.0; ...
%           -1.0,   -2.0,    2.0,   -4.0; ...
%            0.50,  15.0,   -4.0,   21.0];
%       a = [0.1;   0.6;   0.2;   0.4];
%      [Z,L,D,Qz,z] = zt(Qa,a);
%
%    Should return:
%      Qz = [13.0,   -2.0,   -1.5,   0.0; ...
%            -2.0,   10.0,    0.0,   0.0; ...
%            -1.5,    0.0,    4.25,  1.0; ...
%             0.0,    0.0,    1.0,   2.0];
%
%      Z = [0,     0,     1,     0; ...
%           0,     1,     0,     0; ...
%           2,    -1,     1,     1; ...
%           1,    -1,     0,     0];
%
%      z = [0.8;  0.0;  0.3;  0.2];
%
%   See also STUN, EBS, ILS, LTDL.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This function pertains to Section 3.1 in this book.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $

% This routines is based on the decorrelation algorithm
% in the lambda toolbox developed at Delft University of Technology,
% distributed as:
% Function.: decorrel
% Date.....: 19-MAY-1999
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
%
% Changes with respect to that version are mainly related to
% formating, added help, and error checking.



% --- Check input ----------------------------------------------
error(nargchk(1,2,nargin))
if (diff(size(Q)) ~= 0)
  error('Q must be square');
end
if (nargin == 2)
  if ((size(afloat,1)~=size(Q,1)) | (size(afloat,2) ~= 1))
    error('afloat must be a standing vector corresponding to Q');
  end
end



% --- Initialization -------------------------------------------
[L,D] = ltdl(Q);% decomposition Q=L.'*diag(D)*L
n     = size(Q,1);
Zti   = eye(n);
i1    = n - 1;



% --- Guarantee boolean support also for Matlab 5.3 ------------
% --- In later version "true" is a built-in function -----------
% --- Note that I also use '&' not '&&' in order to be ---------
% --- compatible with Matlab 5.3 -------------------------------
true  = logical(1); 
false = logical(0);
% --- The decorrelation procedure ------------------------------
sw    = true;
while (sw==true)
  i  = n;
  sw = false;
  while (( sw==false ) & (i > 1))
    i = i - 1;
    if (i <= i1)
      for j = i+1:n
        mu = round(L(j,i));
        if (mu ~= 0)
          L(j:n,i)   =   L(j:n,i) - mu * L(j:n,j);
          Zti(1:n,j) = Zti(1:n,j) + mu * Zti(1:n,i);
        end
      end
    end
    delta = D(i) + L(i+1,i)^2 * D(i+1);
    if (delta < D(i+1))
      qlambda      = D(i+1) * L(i+1,i) / delta;
      eta          = D(i) / delta;
      D(i)         = eta * D(i+1);
      D(i+1)       = delta;
      %
      % Following 4 lines executed if (i-1 >= 1 && i+1 <= n)
      % Without this if, Matlab generates empty matrices which
      % are handled well by Matlab.
      if (i>1)
        tmp          = L(i+1,1:i-1) - L(i+1,i) .* L(i,1:i-1);
        L(i+1,1:i-1) = qlambda * L(i+1,1:i-1) + eta * L(i,1:i-1);
        L(i,1:i-1)   = tmp;
        L(i+1,i)     = qlambda;
      end
      %
      % Following 3 lines executed if (i+2 <= n)
      % Without this if, Matlab generates empty matrices which
      % are handled well by Matlab.
      if (i+2<=n)
        tmp          = L(i+2:n,i);
        L(i+2:n,i)   = L(i+2:n,i+1);
        L(i+2:n,i+1) = tmp;
      end
      %
      tmp          = Zti(1:n,i);
      Zti(1:n,i)   = Zti(1:n,i+1);
      Zti(1:n,i+1) = tmp;
      %
      % --- Update counters ---
      i1           = i;
      sw           = true;
    end
  end
end



% --- Return Z-transformation matrix after check ----------------------
% --- This is not time consuming, since this function is called only --
% --- once for all PS points, because the design matrix is re-used ----
% --- Note: the back-transform of the fixed solution is obtained as ---
% --- afixed = inv(Z.')*zfixed == Zti*zfixed, i.e., consider ----------
% --- returning "Zti" here to save time later -------------------------
Z  = inv(Zti.');
if (max(max(abs(Z-round(Z))))>1d-6)
  error(['error in Z-transform, deviation from integers: ', ...
          num2str(max(max(abs(Z-round(Z)))))]);
end
Z  = round(Z);



% --- Return decorrelated Q -------------------------------------------
if (nargout >= 4)
  ZtQZ = Z.'*Q*Z;
end



% --- Return the decorrelated ambiguities, if requested ---------------
if ((nargin == 2) & (nargout >= 5))
  zfloat = Z.'*afloat;
end

%%%EOF
