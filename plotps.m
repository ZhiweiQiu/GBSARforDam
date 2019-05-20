function h = plotps(X,Y,Z, ZLIM)
%PLOTPS   plot scattered data triplets.
%   PLOTPS(X,Y,Z) plots the points using the current colormap
%   to the current figure.
%
%   PLOTPS(X,Y,Z, ZLIM) where ZLIM is a two element vectors
%   clips the Z data to the given range.
%
%   H=PLOTPS(...) returns handle H to the figure.
%
%   Example:
%     x = round(rand(100,1)*10000);
%     y = round(rand(100,1)*10000);
%     z = round(randn(100,1)*10);
%     plotps(x,y,z,[-20,20]);
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
% $Revision: 1.3 $  $Date: 2005/12/07 08:23:06 $



% --- Check input -------------------------------------------
error(nargchk(3,4,nargin))
if (length(X)~=length(Y))
  error('length of X and Y must be equal');
end
if (length(X)~=length(Z))
  error('length of X and Z must be equal');
end
if (nargin<4)
  ZLIM = [floor(min(Z)-eps), ceil(max(Z)+eps)];% default
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
clipz_max = max(Z);% plotting range
if (clipz_min==clipz_max)
  ZC   = ones(size(Z));
else
  ZC   = round(1+(NCOL-1-eps).*(Z-clipz_min)./(clipz_max-clipz_min));% [1:64]
end



% --- Plot points with different colors per point -----------
hold off
for i=1:length(X)
  plot(X(i),Y(i),'s', ...
    'Color',cmap(ZC(i),:), ...
    'MarkerFaceColor',cmap(ZC(i),:), ...
    'MarkerSize',4);
  hold on
end
hold off
% --- Add colorbar with corrected annotation ----------------
c = colorbar;
set(c,'YTick',[1,NCOL]);
set(c,'YTIckLabel',[clipz_min,clipz_max]);



% --- Return figure handle if requested ---------------------
if (nargout==1)
  h = gcf;
end

%%% EOF
