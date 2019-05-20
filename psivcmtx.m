function [Q,Fk] = psivcmtx(in,l)
%PSIVCMTX   Variance-covariance matrix model for PSI.
%   Q = PSIVCMTX(K) returns the default vc-matrix for 
%   double-difference observations in a stack of K
%   interferograms using a priori variance factors, i.e.,
%   deg2rad(20)^2 rad^2 for the master SLC image and
%   deg2rad(30)^2 rad^2 for K slave images.
%
%   Q = PSIVCMTX(F) where F is a vector of length K+1
%   returns the vc-matrix using the given factors F.
%
%   QL = PSIVCMTX(K,l) returns the co-factor matrix QL
%   of size (K x K).  l<=K.  Sparse matrix is returned.
%
%   [Q,F] = PSIVCMTX(K) additionally returns the default
%   components of the model in vector F of length K+1 [rad^2].
%
%   The variance component stochastic model is defined as:
%
%          K                  {2*ones(K,K)   if k=0
%     Q = sum F(k+1)*Qk,   Qk=|
%         k=0                 {2*i_k*i_k.'     if k=1:K
%
%   where F(k) are the variance components, Qk are the
%   cofactor matrices, and i_k is vector with a single 1 
%   at position k.
%
%   Examples:
%     Qy = PSIVCMTX(6);% generate a priori vc-matrix
%     Qk = PSIVCMTX(6,0);% returns the first co-factor matrix
%
%   See also STUN, VCE.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This function pertains to Section 2.2.2

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(1,2,nargin))
if (~isreal(in)) 
  error('input must be real');
end
if (prod(size(in))==1) %a priori model
  NIFG = in;
  if (nargin==2)
    % --- Return co-factor matrix Qk ------------------------
    Q = cofactormtx(NIFG,l);
    return
  else
    % --- A priori model ------------------------------------
    Fk    = (30.*pi./180).^2.*ones(NIFG+1,1);
    Fk(1) = (20.*pi./180).^2;
  end
else
  Fk   = in;
  NIFG = length(Fk)-1;% components for SLCs
end



% --- Construct vc-matrix -----------------------------------
Q  = zeros(NIFG);% number of observations of dd time series.
for k=0:NIFG
  Q  = Q + Fk(k+1).*cofactormtx(NIFG,k);
end
%%%EOF



% --- Subfunction to return the cofactor matrix Qk ----------
function Qk = cofactormtx(NIFG,k)
error(nargchk(2,2,nargin))
if (k==0)
  Qk    = 2.*ones(NIFG);% for k=0
else
  i_k    = zeros(NIFG,1);
  i_k(k) = 1;% 0<k<=NIFG
  Qk     = 2.*sparse(i_k*i_k.');% use sparse for speed-up in VCE
end
%%% EOF
