function [Fk,varFk] = psivce(Qy,B,E)
%PSIVCE   Variance Component Estimation for PSI.
%   F = PSIVCE(Q,B,E) returns the estimated variance
%   components of the vc-matrix for K double-difference 
%   observations in a stack of K interferograms.  The
%   square matrix Q of dimension K is the vc-matrix
%   used during a least-squares estimation with design 
%   matrix B [K,NP], with NP the number of estimated 
%   parameters.  Matrix E contains the residuals of the
%   adjusted observations.  If E has more than 1 column,
%   it is assumed that they are uncorrelated and an 
%   variance component estimation is performed for each
%   column, and the returned F are the average.
%
%   [F,VARVAR] = PSIVCE(Q,B,E) additionally returns an
%   estimate for the variance of the estimated (averaged)
%   variance components.
%
%   The variance component stochastic model is defined as:
%
%          K                  {2*ones(K,K)   if k=0
%     Q = sum F(k+1)*Qk,   Qk=|
%         k=0                 {2*i_k*i_k.'     if k=1:K
%
%   where F(k) are the variance components, Qk are the
%   cofactor matrices, and i_k is a vector with a single 1 
%   at position k.
%
%   The variance factors are estimated using
%     F = inv(N)*r;
%   Where
%     r(k+1)     = e.'*inv(Q)*Qk*inv(Q)*e
%     N(k+1,l+1) = trace(inv(Q)*PB*Qk*inv(Q)*PB*Ql)
%   Where
%     e  is a column of E, and
%     PB = eye(K)-B*inv(B.'*inv(Q)*B)*B.'*inv(Q)
%
%   The vc-matrix of the variance factors is
%     D{F} = 2*inv(N)
%   returned is the diagonal of this matrix, if requested.
%   Due to the bad precision of a single estimation of the
%   variance factors it is recommended to average hundreds
%   of independent estimates.  In the STUN algorithm,
%   the estimates are independent because each point is
%   used only once to form a double difference.
%     
%   See also STUN, VCE.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% It pertains to Section 4.3 and Appendix 2.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(3,3,nargin))
NIFG = size(E,1);% number of observations per estimation
NARC = size(E,2);% number of performed estimations
if (size(Qy,1)~=NIFG | size(Qy,2)~=NIFG)
  error('Q must be square, and have same number of lines as e')
end
if (size(B,1)~=NIFG)
  error('B must have same number of lines as Q and e')
end



% --- For each estimation (given least-squares residuals), -----
% --- estimate the components independently in a loop ----------
% --- Estimations should have been performed on arcs -----------
invQy    = inv(Qy);% speed-up; could make it input.
PB       = eye(NIFG)-B*inv(B.'*invQy*B)*B.'*invQy;% orthogonal projector
invQy_PB = invQy*PB;% speed-up
% --- Initialize output matrices -------------------------------
r        = zeros(NIFG+1,1);% number of components is number of SLCs
N        = zeros(NIFG+1,NIFG+1);% number of components
F        = zeros(NIFG+1,1);% sum of estimated components
varF     = zeros(NIFG+1,1);% sum of estimated variance
for arc=1:NARC
  % --- Report current arc, since it seems rather slow ---------
  if (mod(arc-1,20)==0 & NARC>1) 
    disp(['VCE: currently at arc ', num2str(arc), ' of ', num2str(NARC)]);
  end
  % --- Obtain ls residual for this arc and estimate -----------
  e = E(:,arc);
  for k=0:NIFG
    Qk     = psivcmtx(NIFG,k);
    r(k+1) = e.'*invQy*Qk*invQy*e;
    invQy_PB_Qk_invQy_PB = invQy_PB*Qk*invQy_PB;% speed-up
    for l=0:NIFG
      Ql         = psivcmtx(NIFG,l);
      N(k+1,l+1) = trace(invQy_PB_Qk_invQy_PB*Ql);
    end
  end
  invN = inv(N);% better use cholesky for solution
  F    = F    + invN*r;% add this solution
  varF = varF + 2.*diag(invN);
end
Fk = F ./ NARC;% return final estimate as the average.
if (nargout==2)
  varFk = varF ./ NARC ./ (NARC-1);% variance of N independent samples
end

%%% EOF
