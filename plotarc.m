function h = plotarc(X,Y,IDX_from,IDX_to,Z,ZLIM)
%PLOTARC   Plot scattered data points and arcs.
%   PLOTARC(X,Y,IDXFROM,IDXTO) plots the points and arcs indicate
%   by the index vectors to the current figure.
%
%   PLOTARC(X,Y,IDXFROM,IDXTO,Z) where Z is a vector the same
%   length of IDXFROM and IDXTO, uses Z to plot the arcs in
%   color.
%
%   PLOTARC(..., ZLIM) where ZLIM is a two element vectors
%   clips the Z data to the given range.
%
%   H=PLOTARC(...) returns handle H to the figure.
%
%   Example:
%     x    = round(rand(100,1)*10000);
%     y    = round(rand(100,1)*10000);
%     idxf = 1:10;
%     idxt = idxf+1;
%     z    = randn(1,10);
%     plotarc(x,y,idxf,idxt,abs(z),[0 2]);
%
%   See also STUN.

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
% $Revision: 1.2 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(4,6,nargin))
if (nargin<5)
  Z = ones(size(IDX_from));% default same color
end
if (nargin<6)
  ZLIM = [floor(min(Z)-eps), ceil(max(Z)+eps)];% default
end

if (length(X)~=length(Y))
  error('length of X and Y must be equal');
end
if (length(IDX_from)~=length(IDX_to))
  error('length of IDXFROM and IDXTO must be equal');
end
if (length(IDX_from)~=length(Z))
  error('length of IDXFROM and Z must be equal');
end
if (length(ZLIM)~=2)
  error('ZLIM should be a length 2 vector');
end
Z = [Z(:);ZLIM(1);ZLIM(2)]; % add dummies for scaling
Z(find(Z<=ZLIM(1))) = ZLIM(1);
Z(find(Z>=ZLIM(2))) = ZLIM(2);



% --- Create the plot ---------------------------------------
cmap = colormap;% current map
NCOL = size(cmap,1);
% --- Scale Z to range 1:NCOL ---
clipz_min = min(Z);% plotting range
clipz_maz = max(Z);% plotting range
ZC   = round(1+(NCOL-1-eps).*(Z-clipz_min)./(clipz_maz-clipz_min));% [1:64]



% --- Plot arcs with different colors -----------------------
for arc=1:length(IDX_from)
  plot([X(IDX_from(arc)), X(IDX_to(arc))], ...
       [Y(IDX_from(arc)), Y(IDX_to(arc))], ...
      'Color',cmap(ZC(arc),:), ...
      'LineWidth',2);
  hold on
end
plot(X,Y,'ks', 'MarkerFaceColor', [0 0 0], 'MarkerSize',2.0);
hold off



% --- Add colorbar with corrected annotation ----------------
c = colorbar;
set(c,'YTick',[1,NCOL]);
set(c,'YTIckLabel',[clipz_min,clipz_maz]);



% --- Return figure handle if requested ---------------------
if (nargout==1)
  h = gcf;
end

%%% EOF
