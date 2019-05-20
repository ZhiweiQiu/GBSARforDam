%VCEDEMO   Basic demonstration of Variance Component Estimation
%   This script generates unwrapped phase in a stack 
%   of interferograms with known noise levels.  These
%   phase observations are used to estimate DEM error
%   and displacement rate using an a priori stochastic 
%   model.  The least-squares residuals are then
%   used during a VCE to estimate the "true" variance
%   components.
%   Finally, they are used to construct a better vc-matrix
%   which in turn is used to perform a final parameter estimation.
%
%   Note that for this demonstration it is assumed that
%   the unwrapped phase is available.
%   
%   Random data are generated each time this script is 
%   executed.
%
%   Example:
%     vcedemo
%
%   See also STUN, ILSDEMO1D, ILSDEMO1DB.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This script demonstrates the theory in Section 4.3 and Appendix A.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.6 $  $Date: 2005/12/09 16:05:34 $

more off
disp('VCEDEMO   Basic demonstration of Variance Component Estimation');
disp('   This script generates unwrapped phase in a stack');
disp('   of interferograms with known noise levels.  These');
disp('   phase observations are used to estimate DEM error');
disp('   and displacement rate using an a priori stochastic ');
disp('   model.  The least-squares residuals are then');
disp('   used during a VCE to estimate the "true" variance');
disp('   components.');
disp(' ');
disp('   Note that for this demonstration it is assumed that');
disp('   the unwrapped phase is available, and a direct least-');
disp('   squares estimation is performed.  In practice the');
disp('   ILS estimator is used with wrapped data.');
disp(' ');


% --- Remove the reference results if they exist ---------
if (exist('my_Fk')) clear my_*; end


% --- Simulate or load the data --------------------------
disp(' ');
disp('===========================================================');
disp('Enter "1" to load an example to check your installation');
disp('      "-" to generate random input');
disp('===========================================================');
disp(' ');
key = input('Please enter: ','s');
if (isempty(key)) key = '-'; end
switch (lower(key))
  case '1',
    if (~exist('vcedemodat.mat','file'))
      error('   file "vcedemodat.mat" not found');
    end
    load vcedemodat;
    % This data set was created randomly, then saved was:
    % my_Fk     = Fk;
    % my_Qy     = Qy;
    % my_allEST = allEST;
    % my_varFk  = varFk;
    % save vcedemodat NPS NIFG stdNSE Fk_true X Y bperp btemp acq_times PHASE simtopo simdefo phi_ne my_*
    disp('Loading data set "vcedemodat.mat"');
    disp('This is an example to demonstrate a successful estimation.');
    disp('The reference results are stored in the variables:');
    whos my_*
  otherwise,
    rand('state',sum(100*clock));% change state of random number generator
    randn('state',sum(100*clock));% change state of random number generator
    NPS       = 400;
    NIFG      =  30;
    stdNSE    = 20+30*rand(NIFG+1,1);% "true" noiselevels
    stdNSE(1) = 15+10*rand;% set master bit smaller
    Fk_true   = (stdNSE.*pi./180).^2;
    [X,Y]     = simpos(NPS);
    [bperp, btemp, acq_times] = simacq(NIFG);
    [PHASE, simtopo, simdefo, phi_noise] = ...
      simphi(bperp, btemp, X,Y, 50, 20, stdNSE);
end


% --- Plot the simulated data ----------------------------
disp(['Number of PS points:      ', num2str(NPS)]);
disp(['Number of interferograms: ', num2str(NIFG)]);
figure(1)
plotps(X,Y,simtopo);
title('Simulated DEM error');
figure(2)
plotps(X,Y,simdefo);
title('Simulated Deformation Rate');



% --- Perform an LS estimation with standard stoch. model -------
disp('Setting up design matrix for DEM error and lin.defo');
% --- Use same height conversion factor for all points ----
wavelength = 0.056;
slantrange = 850000;
inc_angle  = 23.0*pi/180;% [rad]
KK         = -4*pi/wavelength;
h2p        = KK.*bperp./(slantrange*sin(inc_angle));% [1/m]
v2p        = KK*btemp*1e-3;% [y/mm]
B          = [h2p, v2p];% design matrix

% --- A priori vc-matrix of the dd observations -----------
[Qy, Fk_init] = psivcmtx(NIFG);%

% --- Now form the double-differences at the arcs ---------
% --- Normally, the distances should be similar, but in ---
% --- this case that does not matter, because the noise ---
% --- does not depend on distance. ------------------------
NARC        = NPS/2;% use each point once
disp(['Using ', num2str(NARC), ' independent arcs']);
allEST_init = zeros(2,NARC);
allRES      = zeros(NIFG,NARC);
PROJ_LS     = inv(B.'*inv(Qy)*B)*B.'*inv(Qy);
disp(' ');
disp('Estimate the parameters a priori vc-matrix and unwrapped data.');
disp('  For PSI, actually, the ILS estimator is used, and the wrapped data.');
disp('  This does not affect the principle of Variance Component');
disp('  Estimation, which is demonstrated here.');
disp('  However, note that the standard deviation for uniform noise in the');
disp('  principal interval [-pi,pi) is ~104 degrees, i.e., even if an interferogram');
disp('  contains random noise, the estimated standard deviation using the wrapped');
disp('  residuals will not be larger than 104 degrees!');
disp('Please <Press a key>');
pause
idx_to   = 1:NARC;% PS index
idx_from = NARC+idx_to;% PS index
for arc=1:NARC
  y    = PHASE(:,idx_to(arc))-PHASE(:,idx_from(arc));%dd phase
  xhat = PROJ_LS*y;% estimated DEM error, defo differences
  yhat = B*xhat;
  ehat = y-yhat;
  allRES(:,arc) = ehat;
  allEST_init(:,arc) = xhat;% store for check only
end
disp(' ...done');
disp(' ');

% --- Check estimated parameters with simulated input ----
allDEM_true  = simtopo(idx_to)-simtopo(idx_from);%truth at arcs
allDEFO_true = simdefo(idx_to)-simdefo(idx_from);% truth at arcs
figure(3)
plotps(X(idx_to),Y(idx_to),allEST_init(1,:)-allDEM_true);
title('Error in estimated DEM error at arcs (a priori model)');
figure(4)
plotps(X(idx_to),Y(idx_to),allEST_init(2,:)-allDEFO_true);
title('Error in estimated Displacement Rate at arcs (a priori model)');



% --- Estimate the variance components -------------------
disp('Estimate the variance components (can take some time)');
disp('Please <Press a key>');
pause
[Fk, varFk] = psivce(Qy,B,allRES);
disp(' ...done');
disp(' ');



% --- Compare estimated variance components --------------
xx = 0:NIFG;
figure(5);
plot(xx, 180/pi*sqrt(Fk_true), 'bs', ...
     xx, 180/pi*sqrt(Fk_init), 'y.', ...
     xx, 180/pi*sqrt(Fk),      'r+')
set(gca,'XLim',[-2,NIFG+2], 'YLim',[0,80]);
hold on
errorbar(xx, 180/pi*sqrt(Fk), 180/pi*sqrt(varFk), 'r+');% errors
hold off
title('Variance Components Estimation');
xlabel('SLC number');
ylabel('sqrt(varf) [deg]');
legend('Simulated ("true")', ...
       'A priori', ...
       'Estimated', 1)
disp(' ');
disp('---> Figure 5: Overview of estimated variance factors.');
disp(' ');



% --- Estimate parameters again with new model ------------------
% --- Then second iteration of variance components estimation ---
disp(' ');
disp('Estimate the parameters again, using the a posteriori vc-matrix.');
disp('Please <Press a key>');
pause
Qy      = psivcmtx(Fk);% a posteriori vc-matrix for double-diff obs.
allEST  = zeros(2,NARC);
PROJ_LS = inv(B.'*inv(Qy)*B)*B.'*inv(Qy);
for arc=1:NARC
  y    = PHASE(:,idx_to(arc))-PHASE(:,idx_from(arc));%dd phase
  xhat = PROJ_LS*y;% estimated DEM error, defo differences
  yhat = B*xhat;
  ehat = y-yhat;
  allEST(:,arc) = xhat;% Final parameter estimation
end
disp(' ...done');
disp(' ');

% --- Check estimated parameters with simulated input ----
figure(6)
plotps(X(idx_to),Y(idx_to),allEST(1,:)-allDEM_true);
title('Error in estimated DEM error at arcs (a posteriori model)');
figure(7)
plotps(X(idx_to),Y(idx_to),allEST(2,:)-allDEFO_true);
title('Error in estimated Displacement Rate at arcs (a posteriori model)');
disp(' ');
% --- Compute statistic: with model one is better? --------
disp('Compute statistic: with model one is better?');
disp('--------------------------------------------');
std_dem1  = std(allEST_init(1,:)-allDEM_true);
std_defo1 = std(allEST_init(2,:)-allDEFO_true);
std_dem2  = std(allEST(1,:)-allDEM_true);
std_defo2 = std(allEST(2,:)-allDEFO_true);
disp( '  A priori model/variance factors:');
disp(['    Std.(error estimated topographic term):  ', num2str(std_dem1),  ' [m]']);
disp(['    Std.(error estimated displacement term): ', num2str(std_defo1), ' [mm/y]']);
disp( '  A posteriori model/variance factors:');
disp(['    Std.(error estimated topographic term):  ', num2str(std_dem2),  ' [m]']);
disp(['    Std.(error estimated displacement term): ', num2str(std_defo2), ' [mm/y]']);
disp(' ');
disp('Even if in this demonstration the improvement may be relatively small,');
disp('the VCE is important for ILS, because');
disp('  1. Incorrectly processed interferograms can be detected automatically.');
disp('  2. An over-optimistic assumption on quality of some phase values hamper the');
disp('     non-linear inversion.  Using the correct weights increases the probability');
disp('     of finding the correct integer ambiguities.');
disp('  3. Realistic weighting improves the quality of the estimated parameters.');
disp('  4. Realistic quality description of the estimated parameters is possible.');
disp(' ');



% --- Check with my results if they were loaded from disk --------------
if (exist('my_Fk'))
  if (max(abs(my_Fk-Fk))>1d-3 | ...
      max(abs(my_Qy(:)-Qy(:)))>1d-3 | ...
      max(abs(my_allEST(:)-allEST(:)))>1d-3 | ...
      max(abs(my_varFk-varFk))>1d-3)
    error('difference detected between your results and the reference results!')
  else
    disp('Your results are identical to the reference results!');
  end
end



% --- Say goodbye ------------------------------------------------
disp(' ');
disp('End of demonstration');
disp(' ');
more on
%%%EOF
