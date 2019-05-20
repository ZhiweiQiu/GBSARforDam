function [PHASE, simtopo, simdefo, phi_noise, phi_aps] = ...
  simphi(bperp, btemp, X,Y, maxDEM, stdDEFO, stdNSE, maxAPS, fracdim)
%SIMPHI   Simulate unwrapped interferometric phase in single master stack.
%   PHASE = SIMPHI(BPERP,BTEMP,X,Y) returns a matrix with 
%   interferometric phase values of a single master stack
%   at the positions (X,Y).
%
%   PHASE = SIMPHI(...,maxDEM) simulated the topographic phase
%   using a uniform distribution between [-maxDEM,maxDEM]/2, i.e.,
%   the maximum height difference between two points equals maxDEM.
%   Differential interferograms are simulated; this term
%   corresponds to a DEM error.  A constant height-to-phase
%   conversion factor is assumed.
%   If maxDEM is a vector of the same length as X and Y, then the
%   given values are used directly to simulate the topographic phase.
%
%   PHASE = SIMPHI(...,maxDEM, STD_DEFO) simulates the displacement
%   using a linear model with rates given by STD_DEFO in [mm/y].
%   If STD_DEFO is a vector of the same length as X and Y, then the
%   given values are used directly to simulate the displacement phase.
%
%   PHASE = SIMPHI(...,maxDEM,STD_DEFO, STD_NSE) add noise
%   with a standard deviation STD_NSE [deg] to the unwrapped phase.
%   Noise is simulated on master and slaves.
%   If stdNSE is a vector of length NIFG+1, it is taken as the
%   standard deviation for each SLC image (master is first).
%
%   PHASE = SIMPHI(...,STD_NSE, maxAPS) simulates the 
%   atmospheric phase using a fractal with fractal dimension 
%   2.67 and a maximum variation of maxAPS [rad].  The atmospheric 
%   delays during master and slave acquisitions are simulated and
%   the returned phase is the difference.  If maxAPS is a vector of
%   length NIFG+1, the corresponding value is used for each simulated
%   fractal.
%
%   PHASE = SIMPHI(..., maxAPS, D) uses fractal dimension D to 
%   generate the fractal surfaces.  If D is a vector of length NIFG+1
%   the given fractal dimension is used to generate the individual 
%   fractal surface.
%
%   [PHASE, simtopo, simdefo, phi_aps, phi_nse] = SIMPHI(...) 
%   additionally returns the simulated DEM error, the simulated
%   deformation rates, the atmospheric phase, and the noise
%   phase
%
%   The topographic phase is simulated as:
%     p = -4*pi/wavelength * bperp/(slantrange*sin(inc_angle)) * H;
%   The displacement phase is simulated as:
%     p = -4*pi/wavelength * btemp * 1e-3 * V;
%
%   Typical ERS parameters are used, i.e., wavelength = 0.056 [m],
%   slant-range of 850000 [m], the incidence angle is 23.0 [deg].
%
%   A random bias should be added to the phase in each SLC image to 
%   simulate different offsets due to orbit inaccuracies and in
%   atmospheric delay.  By always considering arcs, the phase
%   differences between points, this bias is always removed 
%   from the double-difference observation, and not simulated here.
%
%   Examples:
%     [X,Y] = simpos(200);
%     [bperp, btemp, acqt] = simacq(30);
%     [PHASE, simtopo, simdefo] = simphi(bperp, btemp, X,Y);
%     figure(1); hist(simtopo);
%     title('histogram of simulated DEM error [m]');
%     figure(2); hist(simdefo);
%     title('histogram of simulated Displacement rates [mm/y]');
%
%   Example to simulate topographic phase with noise at a single point:
%     [PHASE, simtopo, simdefo, phi_noise] = simphi(bperp, btemp, 1,1, 80, 0, 20, 0);
%     figure(3)
%     plot(bperp, PHASE-phi_noise, 'g.', ...
%          bperp, PHASE,'r+', bperp,wrap(PHASE),'b+');
%     title('simulated topographic phase example');
%     xlabel('perpendicular baseline [m]');
%     ylabel('phase [rad]');
%
%   Example to simulate a fractal DEM:
%     [XX,YY] = meshgrid(linspace(1,10000,128));% size [128 x 128]
%     fracdim = 2.2;% for topography on earth
%     [topophase,q,q,q, DEM] = simphi(1, 1, XX(:), YY(:), 0, 0, 0, [1000,0], fracdim);
%     DEM = reshape(DEM,size(XX));
%     DEM = DEM + 200;
%     DEM(find(DEM<0)) = 0;% sea
%     j = jet(64); topomap = [(j(1,:));j(64:-1:15,:)]; colormap(topomap)
%     mesh(DEM); view(-15,60)
%   Note that the output DEM can be used as input during a second run, e.g., 
%   to generate a stack of interferometric phase of spatially correlated 
%   displacement rates, together with temporally changing fractal atmospheres.
%
%   Example to simulate atmospheric signal:
%     [XX,YY] = meshgrid(linspace(1,10000,128));% size [128 x 128]
%     fracdim = 2.67;% atmosphere
%     phi_aps = simphi(1, 1, XX(:), YY(:), 0, 0, 0, [pi,0], fracdim);
%     phi_aps = reshape(phi_aps,size(XX));
%     imagesc(phi_aps); colorbar
%
%   See also STUN, SIMPOS, SIMACQ, ILSDEMO2D.

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
% $Revision: 1.4 $  $Date: 2005/12/07 08:23:06 $

% --- subfunctions in this file ---------------------
% --- function [phi_topo,DEMerror] = simphitopo(bperp,NPS,maxDEM)
% --- function [phi_defo,lindefo]  = simphidefo(btemp,NPS,stdDEFO)
% --- function aps   = simaps(X,Y,NIFG,maxAPS,D)
% --- function noise = simnoise(NIFG,NPS,stdNSE)



% --- Check input -----------------------------------
%         1      2    3 4    5       6        7      8        9
%simphi(bperp, btemp, X,Y, maxDEM, stdDEFO, stdNSE, maxAPS, fracdim)
% ---------------------------------------------------
error(nargchk(4,9,nargin))
if (length(bperp) ~= length(btemp))
  error('bperp and btemp must have same length');
end
if (length(X) ~= length(Y))
  error('X and Y must have same length');
end
if (nargin<5)
  maxDEM = 40;% i.e., between -20:20 [m] of points
end
if (nargin<6)
  stdDEFO = 10.0;% [mm/y] stddev of lin.defo of points
end
if (nargin<7)
  stdNSE = 20.0;% [deg] default noise level on each SLC
end
if (nargin<8)
  maxAPS = 0.0;% [rad] max variation all simulated atmosphere
end
if (nargin<9)
  fracdim = 2.67;% fractal dimension of 2D surfaces
end



% --- Initialize variables --------------------------
NPS  = length(X);
NIFG = length(bperp);



% --- Simulate signal phase -------------------------
% --- Simulate topographic phase --------------------
[phi_topo, simtopo] = simphitopo(bperp,NPS,maxDEM);

% --- Simulate displacement phase -------------------
[phi_defo, simdefo] = simphidefo(btemp,NPS,stdDEFO);

% --- Simulate interferometric atmospheric phase at PS points -
phi_aps = simaps(X,Y,NIFG, maxAPS, fracdim);

% --- Simulate unwrapped interferometric noise at PS points -
phi_noise = simnoise(NIFG,NPS,stdNSE);

% --- Finally, add all components ---------------------------
%PHASE  = phi_topo + phi_defo + phi_aps + phi_noise;%for satellite
PHASE  = phi_defo + phi_aps + phi_noise; %for GBSAR
%%%EOF







% --- Subfunction simphitopo -----------------------------
% --- DEM error uniform error ---
function [phi_topo,DEMerror] = simphitopo(bperp,NPS,maxDEM)
error(nargchk(1,3,nargin))
if (nargin<2)
  NPS = 1;% [m] default
end
if (nargin<3)
  maxDEM = 40.0;% [m] default
end
if (max(maxDEM)<=0)
  phi_topo = 0;
  DEMerror  = zeros(NPS,1);
  return;
end
if (length(maxDEM)~=1 & length(maxDEM)~=NPS)
  error('maxDEM should specify DEM for each PS');
end
% --- Get height per point --------------------------------
if (length(maxDEM)==NPS)
  disp('Using given input DEM errors');
  DEMerror   = reshape(maxDEM,1,NPS);% force lying
else
  rand('state',sum(100*clock));% change state of random number generator
  DEMerror   = maxDEM*rand(1,NPS)-maxDEM./2;
end
% --- Initialize variables --------------------------------
wavelength = 0.0178;
slantrange = 1250;
inc_angle  = 23.0*pi/180;% [rad]
KK         = -4*pi/wavelength;
% --- Use same height conversion factor for all points ----
h2p        = KK.*bperp./(slantrange*sin(inc_angle));
phi_topo   = h2p*DEMerror;% interferometric phase
%EOF



% --- Subfunction simphidefo -----------------------------
% --- lin. displacement rate simulated with normal -------
% --- distribution ---------------------------------------
% --- if nargout ==2 then lindefo, if nargout==3 then...
function [phi_defo,lindefo] = simphidefo(btemp,NPS,stdDEFO)
error(nargchk(1,3,nargin))
if (nargin<2)
  NPS = 1;% [m] default
end
if (nargin<3)
  stdDEFO = 10.0;% [mm/y] default
end
if (max(stdDEFO)<=0)
  phi_defo = 0;
  lindefo  = zeros(NPS,1);
  return;
end
if (length(stdDEFO)~=1 & length(stdDEFO)~=NPS)
  error('stdDEFO should specify DEM for each PS');
end
% --- Get displacement rate per point ---------------
if (length(stdDEFO)==NPS)
  disp('Using given input displacement rates');
  lindefo    = reshape(stdDEFO,1,NPS);% force lying
else
  randn('state',sum(100*clock));% change state of random number generator
  lindefo    = stdDEFO.*randn(1,NPS);% per PS point
end
% --- Initialize variables --------------------------
wavelength = 0.0178;
KK         = -4*pi/wavelength;
v2p        = KK*btemp*1e-3;% [y/mm]
phi_defo   = v2p*lindefo;
% --- more complicated model if more output
if (nargout>=3)
  warning('todo');
end
%EOF



% --- Subfunction simaps -------------------------------
% --- interferometric atmospheric signal at x,y --------
% --- A 2d fractal with fractal dimension D is simulated
% --- related to the exponential beta of a power law as
% --- beta = 7-2D;  for atmosphere beta=5/3 --> D=2.67
% --- using fractal surface with D=5/3 -----------------
% --- Based on code by R.Hanssen -----------------------
function phi_aps = simaps(X,Y,NIFG,maxAPS,D)
error(nargchk(2,5,nargin))
if (nargin<3)
  NIFG = 1;% default single atmosphere
end
if (nargin<4)
  maxAPS = 2;%[rad] default variation [-maxAPS,maxAPS]/2
end
if (nargin<5)
  D = 2.67;% default fractal dimension
end
if (max(maxAPS)<=0)
  phi_aps = 0;
  return;
end
if (length(maxAPS)~=1 & length(maxAPS)~=NIFG+1)
  error('maxAPS should specify noise for each SLC');
end
if (length(maxAPS)==1)
  maxAPS = repmat(maxAPS,NIFG+1,1);% make it a vector
end
if (length(D)~=1 & length(D)~=NIFG+1)
  error('Fractal dimension D should be specified for each SLC');
end
% --- Initialize variables -------------------------
NPS     = length(X);
N       = 512;% represents 10 x 10 km2 area
if (NPS<10000) 
  N = 256;
end
% --- For fractal generation -----------------------
x          = 0.25+(-N/2:N/2-1);% avoid zero distance
[XX,YY]    = meshgrid(x);
beta       = 7-2.*D(1);% 1d-power law coefficient (master)
% --- beta+1 is used as beta, since, the power exponent
% --- is defined for a 1D slice of the 2D spectrum:
beta       = beta+1;
k          = (sqrt(XX.^2 + YY.^2)).^(beta/2);
%
% --- Index where to read fractal ------------------
max_X      = 10000;
max_Y      = 10000;
scX        = 1+round(X*((N-1)/max_X));
scY        = 1+round(Y*((N-1)/max_Y));
idx        = sub2ind([N,N],scY,scX);
phi_aps    = zeros(NIFG,NPS);
master_aps = zeros(1,NPS);
% --- Create the fractals for each acquisition -----
for slc=0:NIFG
  % --- Check whether to adapt maxAPS ---
  CLIM = [-maxAPS(slc+1)/2, maxAPS(slc+1)/2];% output range
  % --- Check whether to adapt beta -----
  if (length(D)==NIFG+1 & slc>0)
    beta = 7-2.*D(slc+1);% 1d-power law coefficient (slave)
    % --- beta+1 is used as beta, since, the power exponent
    % --- is defined for a 1D slice of the 2D spectrum:
    beta = beta+1;
    k    = (sqrt(XX.^2 + YY.^2)).^(beta/2);
  end
  % --- Create this fractals ------------
  rand('state',sum(100*clock));% change state of random number generator
  h = rand(N);% start with random surface
  H = fftshift(fft2(h));
  H = H ./ k;% scale spectrum with power law
  h = abs(ifft2(H));% resulting fractal surface
  h = scale_to(h,CLIM);%
  % --- Get aps at (x,y) ----------------
  this_aps = reshape(h(idx), 1, NPS);% ensure lying vector
  if (slc==0)
    master_aps(:) = this_aps;
  else
    phi_aps(slc,:) = master_aps - this_aps;
  end
end
%EOF



% --- Subfunction scale_to: new limits -------------
function out = scale_to(in,CLIM)
error(nargchk(2,2,nargin))
if (prod(size(CLIM))~=2)
  error('input error scale_to');
end
newmin = CLIM(1);
newmax = CLIM(2);
if (abs(newmin-newmax)<1d-6)
  out = repmat(newmin,size(in));
  return;
end
if (newmin>newmax)
  warning('input error, min>max');
end
minin = min(in(:));
maxin = max(in(:));
if (maxin==minin) 
  warning('cannot scale this');
  out = in;
else
  out = (in-minin).*((newmax-newmin)./(maxin-minin))+newmin;
end
%EOF



% --- Subfunction simnoise -----------------------------
% --- Simulate noise in NIFG interferograms at NPS -----
% --- points.  STDDEV is that of the unwrapped phase ---
% --- in the SLC image, pretending that is possible ----
% --- Noise with a normal distribution is simulated ----
% --- If this is wrapped, it does not have a normal ----
% --- distribution -------------------------------------
% --- If stdNSE is a vector NIFG+1 it is taken as std --
% --- for each SLC image.
% --- Noise is simulated for a master image and --------
% --- the interferometric noise is the difference ------
function phi_noise = simnoise(NIFG,NPS,stdNSE)
error(nargchk(1,3,nargin))
NSLC = NIFG + 1;
if (nargin<2)
  NPS = 1;% default single point
end
if (nargin<3)
  stdNSE = 20.0;% [deg] default noise on SLC
end
if (max(stdNSE)<=0)
  phi_noise = 0;
  return;
end
if (length(stdNSE)~=1 & length(stdNSE)~=NSLC)
  error('stdNSE should specify noise for each SLC');
end
% --- Simulate unwrapped noise in all acquisitions -----
stdNSE = stdNSE.*pi./180;% [rad]
if (length(stdNSE)==1)
  randn('state',sum(100*clock));% change state of random number generator
  p = stdNSE.*randn(NSLC,NPS);%
else
  p = (stdNSE(:)*ones(1,NPS)).*randn(NSLC,NPS);%
end
% --- subtract master (first row) ----------------------
phi_noise = p(2:NSLC,:)-ones(NIFG,1)*p(1,:);
%EOF

%%% EOF

