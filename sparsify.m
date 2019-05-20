function IDX = sparsify(X,Y, W, F, shiftflag)
%SPARSIFY   Sparsification of points.
%   IDX = SPARSIFY(X,Y) returns an index which contains a sparsified data 
%   set, i.e., X(IDX), Y(IDX) is such that only point resides in cell of 
%   size 500 x 500.  IDX is such that Y(IDX) is increasing.
%
%   IDX = SPARSIFY(X,Y, W) uses grid cells with dimensions specified
%   by W=[height,width].  If W is a scalar the bin size is [SY,SY].
%
%   IDX = SPARSIFY(X,Y,W, F) uses F to select the "best" point if a grid
%   cell contains more than a single point.   The returned points have
%   max(F), i.e., if F is the amplitude of the points, then the index of 
%   the points with the largest amplitude in each grid cell are returned.
%
%   IDX = SPARSIFY(..., S) shifts the grid cells with half their dimensions
%   before sparsification.  If S is 1 bins are shifted in Y, if S is 2 they
%   are shifted in X, if S is 3 in both directions.
%
%   The bins are placed over the area as:
%     [   1: SY,   1:SX]  [   1: SY,SX+1:2SX]  ...
%     [SY+1:2SY,   1:SX]  [SY+1:2SY,SX+1:2SX]  ...
%      ...
%
%   Example:
%     [X,Y] = meshgrid(1:1:40);  X=X(:); Y=Y(:); 
%     F     = rand(size(X));
%     IDX   = sparsify(X,Y,[10,5], F);%
%     plot(X,Y,'bo', X(IDX),Y(IDX),'r+');
% 
%   See also STUN.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This function is related to Chapter 4.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $


% --- Check input ----------------------------------------------
error(nargchk(2,5,nargin))
% --- Check dimensions of input --------------------------------
if (~all(size(X)==size(Y)))
  error('X and Y must have same dimensions');
end
% --- Set defaults ---------------------------------------------
if (nargin < 3)
  W = [500,500];% default
end
if (length(W)==1)
  W = [W,W];
end
if (nargin < 4)
  F = ones(size(X));
end
if (nargin < 5)
  shiftflag = 0;% default
end
if (~all(size(X)==size(F)))
  error('X, Y, F must have same dimensions');
end



% --- Bin the data ---------------------------------------------
X         = X-1;
Y         = Y-1;
% --- Shift the grid (same as shift the data) ------------------
switch (shiftflag)
  case 1, 
    Y = Y + W(1)/2;
  case 2, 
    X = X + W(2)/2;
  case 3, 
    Y = Y + W(1)/2;
    X = X + W(2)/2;
end
numbinsX  = 1+floor(max(X)./W(2));
binindex  = floor(X./W(2)) + numbinsX.*floor(Y./W(1));
bins      = unique(binindex);



% --- Search each bin, store index of best_one -----------------
IDX       = zeros(length(bins),1);% allocate output array
for ii=1:length(bins)
  % --- pointer to data present in this bin ---
  idx_this_bin = find(binindex==bins(ii));
  % --- get idx of 'best' one (according to criterium) ---
  [best_one, local_idx_best_one] = max(F(idx_this_bin));
  IDX(ii) = idx_this_bin(local_idx_best_one);
end



% --- Return IDX such that Y(IDX) is increasingly ordered ------
Y2         = Y(IDX);
[q1, IDX2] = sort(Y2);
IDX        = IDX(IDX2);


%%% EOF.

