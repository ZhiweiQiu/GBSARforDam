%STUNDEMO   Basic demonstration of the STUN algorithm.
%   This script simulated the differential interferometric
%   phase at H points in K interferograms.
%
%   1) The precision of the data is estimated.数据精度估计
%   2) A reference network is constructed and the parameters
%      DEM error difference and linear displacement difference
%      are estimated at the arcs. 建立参考网 DEM误差与线形形变估计
%   3) A reference point is selected and the DEM error and 
%      displacement rate are integrated.参考点选取及形变速率整合
%   4) A simplified hypothesis testing is carried out to remove
%      incorrectly estimated parameters.简单假设检验
%
%   See also STUN, VCEDEMO, ILSDEMO1D, ILSDEMO1DB, SIMPHI.

% This function is part of the STUN toolbox which
% accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.
%
% This function pertains to Chapter 4 of this book.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.5 $  $Date: 2005/12/07 08:23:06 $



% --- Control parameters for skipping steps for debugging --------------



% --- Start the demonstration ------------------------------------------
more off
disp('STUNDEMO   Basic demonstration of the STUN algorithm.');
disp(' ');
disp('===========================================================');
disp('Enter "1" to load data set with fractal deformation');
disp('      "2" to load data set with fractal deformation/atmospheres');
disp('      "l" to load a saved data set of a previous run');
disp('      "p" to profile this run');
disp('      "-" to simulate randomly');
disp('===========================================================');
disp(' ');
key = input('Please enter: ','s');
if (isempty(key)) key = '-'; end;
do_profile = 0;
switch (lower(key))
  case '1',
    disp('Loading file "stundemodat1"');
    disp('Data simulated using a fractal displacement field');
    disp('No atmosphere');
    if (~exist('stundemodat1.mat','file'))
      error('   file "stundemodat1.mat" not found');
    end
    load stundemodat1
  case '2',
    disp('Loading file "stundemodat2"');
    disp('Data simulated using a fractal displacement field');
    disp('Same data as "1", with additional atmosphere');
    if (~exist('stundemodat2.mat','file'))
      error('   file "stundemodat2.mat" not found');
    end
    load stundemodat2
  case 'l',
    [FILENAME, PATHNAME] = uigetfile('*.mat', 'load stundemo file');
    file = [PATHNAME, FILENAME];
    if (~exist(file,'file'))
      error(['file ', file, ' does not exist']);
    end
    load(file);
  otherwise,
    % --- Check if profiling is requested ------------------------------
    if (lower(key)=='p')
      do_profile = 1;
      profile clear;
      profile on;
    end
    % --- Simulate the data --------------------------------------------
    rand('state',sum(100*clock));% change state of random number generator
    randn('state',sum(100*clock));% change state of random number generator
    NIFG   = input('Enter number of IFG [30]:          ');
    NPS    = input('Enter number of PS  [1000]:        ');
    qNSE   = input('Enter avg. noise level at SLC [20] deg: ');
    maxAPS = input('Enter max. atmospheric signal [0]  rad: ');
    if (isempty(NIFG))   NIFG=30;  end
    if (isempty(NPS))    NPS=1000; end
    if (isempty(qNSE))   qNSE=20;  end
    if (isempty(maxAPS)) maxAPS=0; end
    % --- Noise per SLC -------------------------------------------------
    if (qNSE<eps) 
      stdNSE = 0;
    else
      stdNSE    = qNSE-10+20*rand(NIFG+1,1);% "true" noiselevels
      stdNSE(1) = 15+10*rand;% set master bit smaller
    end
    % --- Simulate all data ---------------------------------------------
    Fk_true   = (stdNSE.*pi./180).^2;% "true" variance components
    maxDEM    = 50;% [-25,25] meter DEM error, uniformly distributed
    stdDEFO   = 10;% normally distributed, spatially random, displacement rates
    [X,Y]     = simpos(NPS);
    [Bperp, Btemp, acq_times] = simacq(NIFG);
    [PHASE_UW, simtopo, simdefo, phi_noise, phi_aps] = ...
      simphi(Bperp, Btemp, X, Y, maxDEM, stdDEFO, stdNSE, maxAPS);
    PHASE = wrap(PHASE_UW);% observations are wrapped
    disp('Wrapped phase observation simulated');
end


%    是否进行方差分量估计
% --- Ask user if VCE is desired --------------------------------------
NO_VCE = input('Enter "1" to *not* perform VCE [0]:     ');
if (isempty(NO_VCE)) NO_VCE=0; end
if (NO_VCE==1)
  warning('only for quicker testing')
end



% --- Report overview -------------------------------------------------
disp(' ');
disp('Overview of parameters:');
disp('------------------------------');
disp(['  Number of points:                       ', num2str(NPS)]);
disp(['  Number of interferograms:               ', num2str(NIFG)]);
disp(['  Width of area [km]:                     ', ...
  num2str((max(X)-min(X))/1000)]);
disp(['  Height of area [km]:                    ', ...
  num2str((max(Y)-min(Y))/1000)]);
disp(['  Span of Perpendicular baselines [m]:    ', ...
  num2str(max(Bperp)-min(Bperp))]);
disp(['  Span of Temporal baselines [y]:         ', ...
  num2str(max(Btemp)-min(Btemp))]);
disp(['  Max. DEM error [m]:                     ', num2str(max(simtopo))]);
disp(['  Min. DEM error [m]:                    ', num2str(min(simtopo))]);
disp(['  Max. Displacement rate [mm/y]:          ', num2str(max(simdefo))]);
disp(['  Min. Displacement rate [mm/y]:         ', num2str(min(simdefo))]);
disp(['  Average noise level [deg]:              ', num2str(mean(stdNSE))]);


%   对目标点进行稀疏化
% --- Sparsify the points to obtain reference network points --
NO_SPARSIFICATION=0;
if (NO_SPARSIFICATION==1)
  warning('Not performing sparsification by user request')
  IDX_REF = 1:length(X);
else
  F       = rand(size(X));% e.g., inverse amplitude dispersion index
  IDX1    = sparsify(X,Y,[500,500], F,0);%
  IDX2    = sparsify(X(IDX1),Y(IDX1),[500,500], F(IDX1),1);%
  IDX3    = sparsify(X(IDX1(IDX2)),Y(IDX1(IDX2)),[500,500], F(IDX1(IDX2)),2);%
  IDX_REF = IDX1(IDX2(IDX3));% index in X,Y
end
X_REF    = X(IDX_REF);
Y_REF    = Y(IDX_REF);


%     建立狄洛尼三角网及建立弧向量指数
% --- Create a network using Delaunay --------------------------
% --- Create arcs index vectors that specify this network ------
TRI      = delaunay(X_REF,Y_REF);
all_from = [TRI(:,1); TRI(:,2); TRI(:,3)].';
all_to   = [TRI(:,2); TRI(:,3); TRI(:,1)].';
% --- Make a list from-to for all arcs, sort it ---
% --- Make sure from is always less than to -------
qfrom    = min(all_from,all_to);% index in X_REF
qto      = max(all_from,all_to);% index in X_REF
ft       = [qfrom.', qto.'];
ft       = unique(ft,'rows');% remove double arcs
IDX_from = ft(:,1);% index in ?_REF (sorted)
IDX_to   = ft(:,2);% index in ?_REF
clear ft qfrom qto TRI;


%   去除过长的弧线
% --- Remove long arcs -----------------------------
all_arclengths     = sqrt((X_REF(IDX_from)-X_REF(IDX_to)).^2 + ...
                          (Y_REF(IDX_from)-Y_REF(IDX_to)).^2);
max_length         = 5000;
too_long           = find(all_arclengths > max_length);
IDX_from(too_long) = [];% remove these arcs
IDX_to(too_long)   = [];%  +assume all points remain
all_arclengths     = sqrt((X_REF(IDX_from)-X_REF(IDX_to)).^2 + ...
                          (Y_REF(IDX_from)-Y_REF(IDX_to)).^2);

%    记录弧线上为“真值”参数
% --- Store the "true" parameters at the arcs ------
% --- of the reference network ---------------------
allDEM_arc_true   = simtopo(IDX_REF(IDX_to))-simtopo(IDX_REF(IDX_from));
allDEFO_arc_true  = simdefo(IDX_REF(IDX_to))-simdefo(IDX_REF(IDX_from));



% --- Report stats ---------------------------------
NARC = length(IDX_from);
NREF = length(unique([IDX_from;IDX_to]));%
disp(' ');
disp('Overview of reference network:');
disp('------------------------------');
disp(['  Number of points in network:            ', num2str(NREF)]);
disp(['  Number of arcs:                         ', num2str(NARC)]);
disp(['  Average number of arcs per point:       ', num2str(2*NARC/NREF)]);% 1 arc is 2ps
disp(['  Minimum arc length:                     ', num2str(min(all_arclengths))]);
disp(['  Maximum arc length:                     ', num2str(max(all_arclengths))]);
disp(['  Mean arc length:                        ', num2str(mean(all_arclengths))]);
disp(['  Standard deviation arc length:          ', num2str(std(all_arclengths))]);
disp(' ');
disp('Please <Press a key> to continue');
if (do_profile~=1) pause; end;


%  为整合产生一个设计矩阵
% --- Create a design matrix for "integration" ------
% --- do not reduce it for reference point yet ------
LS_DESIGN_C = zeros(NARC, NREF);% chapter 4.4.2
for arc=1:NARC
  LS_DESIGN_C(arc,IDX_from(arc)) = -1; 
  LS_DESIGN_C(arc,IDX_to(arc))   =  1; 
end
NARC_per_point = sum(abs(LS_DESIGN_C));


%  为了进行方差分量估计，从弧线中寻找一些没共用的单点
% --- For VCE, find a set of arcs such that each point is used only once ---
disp(' ');
disp('Obtaining a set of arcs for VCE that do not use common points');
tmp = LS_DESIGN_C;
for p=1:NREF
  used = find(tmp(:,p));% e.g., in rows [1,2]
  tmp(used(2:length(used)),:)=[];% remove all rows that also have this point
end
NARC_vce     = size(tmp,1);
IDX_from_vce = zeros(1,NARC_vce);
IDX_to_vce   = zeros(1,NARC_vce);
for arc=1:NARC_vce
  IDX_from_vce(arc) = find(tmp(arc,:)==-1);
  IDX_to_vce(arc)   = find(tmp(arc,:)==1);
end
clear tmp;% not needed anymore


%   计算每个弧线的距离
% --- Compute distances of each arc for VCE ------------------------
arclengths = sqrt((X_REF(IDX_from_vce)-X_REF(IDX_to_vce)).^2 + ...
                  (Y_REF(IDX_from_vce)-Y_REF(IDX_to_vce)).^2);
disp(['  Number of arcs for VCE:                 ', num2str(NARC_vce)]);
disp(['  Mean arc length for VCE:                ', num2str(mean(arclengths))]);
disp(['  Standard deviation arc length for VCE:  ', num2str(std(arclengths))]);


%    绘制获取影像分布网
% --- Plot acquisitions/network ------------------------------------
figure(1)
  %%% Panel 1: baseline distribution 基线分布
  subplot(2,2,1)
    plot([0;Bperp],acq_times,'r+');% add master
    title('Baseline distribution');
    xlabel('perpendicular baseline');
    ylabel('acquisition time');
  %%% Panel 2: simulated deformation 模拟形变
  subplot(2,2,2)
    plotps(X_REF, Y_REF, simdefo(IDX_REF), 2*[-stdDEFO,stdDEFO]);
    title('Simulated displacement rates');
    xlabel('Ground-range [m]');
    ylabel('Azimuth [m]');
  %%% Panel 3: spatial distribution of points, network 点位空间分布网
  subplot(2,2,3)
    plotarc(X_REF, Y_REF, IDX_from, IDX_to);
    hold on
    plot(X, Y, 'k.', 'MarkerSize',2);
    hold off
    title('Reference Network');
    xlabel('Ground-range [m]');
    ylabel('Azimuth [m]');
  %%% Panel 4: arcs used for VCE  绘制弧段
  subplot(2,2,4)
    plotarc(X_REF, Y_REF, IDX_from_vce, IDX_to_vce);
    title('Arcs for Variance Component Estimation');
    xlabel('Ground-range [m]');
    ylabel('Azimuth [m]');
disp(' ');
disp('Figure 1 show the simulated data');

% --- Plot some atmospheres -------绘制大气相位--------
if (maxAPS~=0)
figure(11)
  %%% Panel 1: baseline distribution 绘制基线分布
  subplot(2,2,1)
    plotps(X, Y, phi_aps(1,:), 0.5*[-maxAPS,maxAPS]);
    title('Simulated atmosphere 1');
  %%% Panel 2: simulated deformation 形变模拟
  if (NIFG>=2)
    subplot(2,2,2)
      plotps(X, Y, phi_aps(2,:), 0.5*[-maxAPS,maxAPS]);
      title('Simulated atmosphere 2');
  end
  %%% Panel 3: spatial distribution of points, network 空间点位分布网
  if (NIFG>=3)
    subplot(2,2,3)
      plotps(X, Y, phi_aps(3,:), 0.5*[-maxAPS,maxAPS]);
      title('Simulated atmosphere 3');
  end
  %%% Panel 4: arcs used for VCE 绘制弧段
  if (NIFG>=4)
    subplot(2,2,4)
      plotps(X, Y, phi_aps(4,:), 0.5*[-maxAPS,maxAPS]);
      title('Simulated atmosphere 4');
  end
disp(' ');
disp('Figure 11 show the simulated atmosphere for IFG=1:4');
end

disp('Please <Press a key> to continue');
if (do_profile~=1) pause; end;



%   为最小二乘整数修正（ILS）建立函数模型
% ------------------------------------------------------------------
% --- Set up functional model for ILS ------------------------------
% ------------------------------------------------------------------
disp('Setting up design matrix for DEM error and lin.defo');
% --- Use same height conversion factor for all points ----
wavelength = 0.056;% [m]
slantrange = 850000;% [m]
inc_angle  = 23.0*pi/180;% [rad]
KK         = -4*pi/wavelength;
h2p        = KK.*Bperp./(slantrange*sin(inc_angle));% [1/m]
v2p        = KK*Btemp*1e-3;% [y/mm]
B          = [h2p, v2p];% design matrix for float parameters
NPM        = size(B,2);% number of float parameters (and pseudo-obs)


%    为ILS建立随机模型
% --- Set up stochastic model for ILS ------------------------------
%    一个双差观测的先验方差分量矩阵构建
% --- A priori vc-matrix of the dd observations --------------------
disp('Qy: set to a priori vc-matrix of the DD-observations');
disp('See section 4.3 of the book');
[Qy, Fk_init] = psivcmtx(NIFG);%


%     方程组规则化并加入观测值残差和方差
% --- Regularization of the system of equations --------------------
% --- Add pseudo-observations and their variances ------------------
disp(' ');
disp('Setting up the model of observation equations:');
disp('   y1 = A1*a + B1*b + e1;  Qy1=D{e1},              (Eq. 1)');
disp('where:');
disp('  y1 [N x 1] vector with N wrapped phase observations.');
disp('  A1 [N x N] design matrix for integer parameters.');
disp('  a  [N x 1] vector of N unknown integer parameters (the ambiguities).');
disp('  B1 [N x 1] design matrix for float parameters.');
disp('  b  [2 x 1] vector of unknown float parameters.');
disp('  e1 [N x 1] vector of measurement noise.');
disp(' Qy1 [N x N] vc-matrix of the noise (D{e} is dispersion of e).');
disp(' ');
disp('This system of equations (Eq. 1) is under-determined,');
disp('i.e., there are more unknown parameters than observations.');
disp('Therefor, (Eq. 1) is regularized using 2 pseudo-observations');
disp('for the float parameters, with value 0 and certain standard deviation.');
disp('  [y1] = [A1] * a + [B1] * b + [e1];  Qy=[D{e1}  0 ],  (Eq. 2)');
disp('  [y2]   [A2]       [B2]       [e2]      [ 0  D{e2}]');
disp(' ');
%    加入伪观测量
% --- Add pseudo-observation -------------------------------
%y1     = y;
%y2     = 0;% zero pseudo-observation
y2      = zeros(NPM,1);% zero pseudo-observation加入伪观测量
%y      = [y1; y2];

%     为整数参数建立设计矩阵
% --- Set up design matrix for integer parameters ----------
A1      = -2*pi*eye(NIFG);
A2      = zeros(NPM,NIFG);
A       = [A1; A2];

%     为浮点参考建立设计矩阵
% --- Set up design matrix for float parameters ------------
B1      = B;%
B2      = eye(NPM);
B       = [B1; B2];


%     建立方差分量矩阵
% --- Set up vc-matrix -------------------------------------
var_pseudo_obs = [0.5*maxDEM, 2*stdDEFO].^2;% set a priori "soft-bounds"建立先验软边界
Qy1     = Qy;
Qy2     = diag(var_pseudo_obs);
Qy      = [Qy1,     zeros(NIFG,NPM); ...
           zeros(NPM,NIFG),    Qy2];


%     利用LAMBDA方法进行估计 快速进行方程组去相关并进行自动估计
% --- Perform the estimation using the LAMBDA method -------
% --- Decorrelation of system of equations for faster and --
% --- more robust estimation -------------------------------
% --- See chapter 3 ----------------------------------------
P_B     = eye(NIFG+NPM) - B * inv(B.'*inv(Qy)*B) * B.'*inv(Qy);
A_      = P_B*A;% reduced design matrix
Qahat   = inv(A_.'*inv(Qy)*A_);%


%     利用Z转换进行去相关 
% --- Decorrelation using Z-transform, function "zt" -------
% --- Speed-up, compute Z-transform outside loop (independent of observations)
[Z,L,D] = zt(Qahat);
Linv    = inv(L);% pre-computed outside of loop
Dinv    = 1./D;% pre-computed outside of loop
invZt   = inv(Z.');% for fixed solution a_check
invLtDL = inv(L.'*diag(D)*L);% pre-computed outside of loop for "ebs"


%     利用解缠数据进行浮点参数估计
% --- 7: Estimate float parameter using unwrapped data -----
% --- Least-squares solution of E{y}=Ax; D{y}=Qy -----------
% --- is given by xhat=inv(A.'*inv(Qy)*A)*A.'*inv(Qy)*y ----
PROJ_LS     = inv(B1.'*inv(Qy1)*B1)*B1.'*inv(Qy1);% for speed


%    利用先验随机模型进行弧段的方差分量估计
% --- Wait for user ----------------------------------------
% --- Then perform a VCE -------------------------------------------
% --- Use these arcs to estimate using a priori stochastic model ---
disp('Set up of matrices completed');
disp(' ');
disp(' ');
disp('Start of estimate the variance components.');
disp(['  First estimate parameters at ', num2str(NARC_vce), ' independent arcs']);
disp('  Using the Extended Bootstrap using the a priori vc-matrix (wrapped data).');
disp('Please <Press a key> to continue');
if (do_profile~=1) pause; end;


%     由弧段双差中估计
% --- Form the double-differences at the arcs ------------------------
allEST_init = zeros(NPM, NARC_vce);
allsqnorm   = zeros(1,NARC_vce);
allRES      = zeros(NIFG,NARC_vce);
for arc=1:NARC_vce
  %        弧段双差相位观测值计算
  % --- 0: Double-difference phase observations at arc ---------------
  y1   = wrap(PHASE(:,IDX_to_vce(arc))-PHASE(:,IDX_from_vce(arc)));
  %y    = [y1;y2];% add pseudo-observations
  %        浮点模糊度计算
  % --- 1: float solution for ambiguities ----------------------------
  afloat = y1./(-2.*pi);% float solution (special case for PSI)
  zfloat = Z.' * afloat;% decorrelate float solution for faster ILS search
  %        为最小二乘整数估计获取边界
  % --- 2: Obtain bound for ILS search, function "ebs" ---------------
  [zfixed,sqnorm] = ebs(zfloat,L,D,1,invLtDL);% 1 candidate, reuse invLtDL
  %        a 修正方法，Z转换的逆
  % --- 4: Fixed solution a_check, Inverse Z-transform ---------------
  afixed = invZt*zfixed;
  %        利用估计模糊度进行数据解缠
  % --- 5: Unwrap data using the estimated ambiguities ---------------
  y_uw   = y1 + 2.*pi.*afixed;
  %        利用解缠数据进行常规的带权最小二乘
  % --- 6: Ordinary weighted least-squares using unwrapped data ------
  bhat   = PROJ_LS*y_uw;% estimated DEM error, defo differences
  yhat   = B1*bhat;
  ehat   = wrap(y_uw-yhat);% residuals are always wrapped, closest solution
  allRES(:,arc)      = ehat;
  allsqnorm(arc)     = sqnorm;
  allEST_init(:,arc) = bhat;% store to get feeling for parameter values
end
disp(' ...done');
disp(' ');

%在这些弧段上进行内估计，估计方差分量并利用他们进行参考网的方差参数估计
%
% --- After the initial estimation at these arcs, estimate the -------
% --- variance components and use them to estimate the parameters ----
% --- at all arcs of the reference network ---------------------------
% --- Estimate the variance components -------------------------------
disp('Estimate the variance components (can take some time)');
disp('Please <Press a key>');
if (do_profile~=1) pause; end;
IDX_ok      = find(allsqnorm<mean(allsqnorm)+std(allsqnorm));
disp(['Using ', num2str(length(IDX_ok)), ' arcs of ', num2str(NARC_vce), ' for estimation of Fk']);
if (NO_VCE==1)
  warning('Not performing VCE by user request');
  Fk = Fk_true;
else
  [Fk, varFk] = psivce(Qy1,B1,allRES(:,IDX_ok));
end
disp(' ...done');


%    核对方差分量是否估计正确
% --- Check if Variance Components were estimated correctly ----------
% --- Since we know their "true" value here --------------------------
errFk   = max(abs(rad2deg(sqrt(Fk_true)-sqrt(Fk))));
errFk_  = mean((rad2deg(sqrt(Fk_true)-sqrt(Fk))));
errFk__ = std(abs(rad2deg(sqrt(Fk_true)-sqrt(Fk))));
disp(['Max. error of estimated standard deviation of data [deg]: ', num2str(errFk)]);
disp(['Mean error of estimated standard deviation of data [deg]: ', num2str(errFk_)]);
disp(['Std. error of estimated standard deviation of data [deg]: ', num2str(errFk__)]);
if (maxAPS~=0 & errFk_<0)
  disp('  A negative mean indicates that the error is estimated larger');
  disp('  than the random noise level, e.g., atmospheric signal is');
  disp('  also accounted for in the estimated variance components.');
end
disp(' ');


%     随机模型的最小二乘整数参数估计
% --- ILS Estimate parameters with this stochastic model --------------
disp(' ');
disp('Estimate the parameters at arcs of network using the a posteriori vc-matrix.');
disp('Please <Press a key>');
if (do_profile~=1) pause; end;


%     更新方差分量矩阵，利用主要参数设定参数边界
% --- Update vc-matrix -------------------------------------------------
% --- Set bounds on parameters using preliminary estimates -------------
% --- Qy1: a posteriori vc-matrix for double-diff observations ---先验的双差观测方差分量矩阵------
if (NO_VCE~=1)
  Qy1     = psivcmtx(Fk);%<-- evaluate estimated components!
end

% --- Qy2: vc-matrix for pseudo-observations --------残差观测方差分量矩阵-------------------
var_pseudo_obs = (1.5*std(allEST_init(:,IDX_ok).')).^2;%
disp(['Std. of pseudo-observations for DEM error, displ. rate: ', num2str(sqrt(var_pseudo_obs))]);
Qy2     = diag(var_pseudo_obs);%<-- derived from preliminary estimates
Qy      = [Qy1,     zeros(NIFG,NPM); ...
           zeros(NPM,NIFG),    Qy2];


%     新模型更新
% --- Update with new model --------------------------------------------
Qahat   = inv(A_.'*inv(Qy)*A_);% functional model same, stochastic model updated


%     利用Z变换进行去相关处理
% --- Decorrelation using Z-transform, function "zt" ---------------------
%     加速：计算跳出循环（独立观测量）
% --- speed-up: computations outside loop (independent of observations) --
[Z,L,D] = zt(Qahat);
Linv    = inv(L);% pre-computed outside of loop 跳出运算预计算
Dinv    = 1./D;% pre-computed outside of loop 
invZt   = inv(Z.');% for fixed solution a_check
invLtDL = inv(L.'*diag(D)*L);% pre-computed outside of loop for "ebs"


%     利用解缠数据估计浮点参数
% --- Projector to estimate float parameter using unwrapped data -------
Qb_hat     = inv(B1.'*inv(Qy1)*B1);% vc-matrix of estimated parameters 
PROJ_LS    = Qb_hat*B1.'*inv(Qy1);% for speed outside loop 快速跳出
disp(['Estimation at ', num2str(NARC), ' arcs of the network']);
disp('Qb_hat, vc-matrix of estimated parameters follows: ');
Qb_hat = Qb_hat;
disp(' ');


%     由弧段双差进行估计
% --- Form the double-differences at the arcs ------------------------
allESTarc   = zeros(NPM, NARC);
allsqnorm   = zeros(1,NARC);
allRESarc   = zeros(NIFG,NARC);
for arc=1:NARC
  if (mod(arc-1,50)==0)
    disp(['  ILS: arc ', num2str(arc), ' of ', num2str(NARC)]);
  end
  %      弧段双差观测相位
  % --- 0: Double-difference phase observations at arc ---------------
  y1   = wrap(PHASE(:,IDX_REF(IDX_to(arc))) - ...
              PHASE(:,IDX_REF(IDX_from(arc))));
  %      加入残差观测
  %y    = [y1;y2];% add pseudo-observations
  % --- 1: float solution for ambiguities ----------------------------
  afloat = y1./(-2.*pi);% float solution (special case for PSI)
  zfloat = Z.' * afloat;% decorrelate float solution (for faster search)
  % --- 2: Obtain bound for ILS search, function "ebs" ---------------
  [zfixed_bs, sqnorm_bs]   = ebs(zfloat,L,D,1,invLtDL);% 1 candidate
  % --- 3: Fixed solution for decorrelated ambiguities using "ils"----
  [zfixed_ils, sqnorm_ils, ierr] = ils(zfloat,Linv,Dinv,sqnorm_bs,1);
  % --- Check ierr status flag ---------------------------------------
  if (bitand(ierr,1))
    warning('not enough candiates found in "ils"');
  end
  if (bitand(ierr,2))
    warning('maxloop reached in "ils", using bootstrapped estimate');
    sqnorm_ils = sqnorm_bs;
    zfixed_ils = zfixed_bs;
  end
  % --- 4: Fixed solution a_check, Inverse Z-transform ---------------
  afixed = invZt*zfixed_ils;
  % --- 5: Unwrap data using the estimated ambiguities ---------------
  y_uw   = y1 + 2.*pi.*afixed;
  % --- 6: Ordinary weighted least-squares using unwrapped data ------
  bhat   = PROJ_LS*y_uw;% estimated DEM error, defo differences
  yhat   = B1*bhat;
  ehat   = wrap(y_uw-yhat);% residuals are always wrapped, closest solution
  allRESarc(:,arc) = ehat;
  allsqnorm(arc)   = sqnorm_ils;
  allESTarc(:,arc) = bhat;% store estimated parameters at arcs
end
disp(' ...done');
disp(' ');


%     真值检核
% --- Check with "true" ------------------------------------------------
% error of estimated parameters at arcs of network
errEST = [allDEM_arc_true;allDEFO_arc_true]-allESTarc;
stdERR = std(errEST.');
disp('Comparing estimated differences with simulated "true" values');
disp(['  Std. of difference with "true" DEM error [m] at arcs: ', ...
  num2str(stdERR(1))]);
disp(['  Std. of difference with "true" Displ.rate [mm/y] at arcs: ', ...
  num2str(stdERR(2))]);
disp(' ');


%     对一些参考点进行简单最小二乘整合
% --- Simple ls-integration with respect to some reference point -------
% --- yields an idea of the parameters at the points. ------------------
% --- Errors with obs, remove lots of arcs inmediately -----------------
% --- Thorough alternative hypothesis testing (DIA) is not performed ---
disp('Least-Squares integration of estimated differences at arcs.');
disp('  If all estimations at the arcs are correct, then the adjustes');
disp('  residuals at the arcs should all be zero, and thus, integration');
disp('  would be path independent.');
IDX_refpnt      = 1;
C               = LS_DESIGN_C;% integration of estimates at arcs
C(:,IDX_refpnt) = [];% remove singularity, DEM_refpnt==0; DEFO_refpnt==0
disp(['Selected reference point with IDX = ', num2str(IDX_refpnt)]);
% --- remove estimated with large sqnorm -------------------------------
% --- This is not a good way of hypothesis testing ---------------------
% --- and no test for points -------------------------------------------
IDX_ok    = find(allsqnorm<mean(allsqnorm)+2*std(allsqnorm)+eps);
y         = allESTarc(:,IDX_ok).';% remove bad arcs
C         = C(IDX_ok,:);% remove bad arcs
IDX_from  = IDX_from(IDX_ok);% reduce for removed arcs
IDX_to    = IDX_to(IDX_ok);% reduce for removed arcs
disp(['Removing ', num2str(NARC-length(IDX_ok)), ' arcs of ', ...
  num2str(NARC),' that had a bad fit before inversion']);
disp('  if many arcs are removed like this, it may lead');
disp('  to points that are not connected with any arc (which should be taken out),')
disp('  or to isolated networks that require additional reference points')
disp('  Here, this is not tested, i.e., the following system of equations');
disp('  may be singular, which will then lead to a crash.');
points_with_no_arcs = find(sum(abs(C))==0);
disp(['  Number of unconnected points:         ', ...
  num2str(length(points_with_no_arcs))]);


%     转置
% --- Perform the inversion --------------------------------------------
%     假设参考错误
% --- Assume that the parameters are uncorrelated ----------------------
Qy_est    = speye(size(C,1));% inv. weights for LS-integration
LS_INT    = inv(C.'*inv(Qy_est)*C)*C.'*inv(Qy_est);
PM_at_REF = LS_INT * y;% xhat


%     计算改正值
% --- Compute the adjusted variables -----------------------------------
yhat      = C*PM_at_REF;% adjusted observations
ehat      = y-yhat;% adjusted residuals at arcs


%     假设检验去除更多的弧段
% --- One should now perform hypothesis testing to remove more arcs ----
% --- However, we continue with the algorithm --------------------------
% --- for demonstration purposes ---------------------------------------
sum_sq_err     = sum(ehat.^2,2);
[q, IDX_worst] = max(sum_sq_err);
disp(' ');
disp('Arc with largest adjusted residuals');
disp(['  refnet point IDX ', ...
  num2str(IDX_from(IDX_worst)), ' ---> to ' num2str(IDX_to(IDX_worst))]);
disp(['  adjusted residual DEM error:    ', num2str(ehat(IDX_worst,1))]);
disp(['  adjusted residual Displ. rate:  ', num2str(ehat(IDX_worst,2))]);
if (sum_sq_err>2)
  disp('Large residual');
  disp('This arc should be removed to obtain a consistent network');
  disp('However, this program will continue');
end


%     再次回来加入参考点
% --- Add the reference point back in again ----------------------------
% --- after this, the order of the points is as in X_REF ---------------
PM_at_REF = [PM_at_REF(1:IDX_refpnt-1,:); ...
             [0,0]; ...
             PM_at_REF(IDX_refpnt:size(PM_at_REF,1),:)];


%     绘制弧段的改正残差值
% --- Plot the adjusted residuals at the arcs --------------------------
figure(2)
  subplot(2,2,1)
    plotarc(X_REF, Y_REF, IDX_from, IDX_to, ehat(:,1),[-10,10]);
    title('Adjusted residuals at arcs for DEM error')
  subplot(2,2,2)
    plotarc(X_REF, Y_REF, IDX_from, IDX_to, ehat(:,2),[-10,10]);
    title('Adjusted residuals at arcs for displ. rate')
  subplot(2,2,3)
    plotps(X_REF, Y_REF, PM_at_REF(:,1), [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Estimated DEM error [m] at PS points')
  subplot(2,2,4)
    plotps(X_REF, Y_REF, PM_at_REF(:,2), [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Estimated Displ. rates [mm/y] at PS points')


%     检验是否获取正确参考
% --- Check if the correct parameters are obtained ---------------------
% --- Since we know their "true" value here ----------------------------
PM_at_REF_true  = [(simtopo(IDX_REF)-simtopo(IDX_REF(IDX_refpnt))).', ...
                   (simdefo(IDX_REF)-simdefo(IDX_REF(IDX_refpnt))).'];
PM_at_REF_error = PM_at_REF_true-PM_at_REF;
disp(['Std. of difference with simulated input: ', num2str(std(PM_at_REF_error))]);
disp(' ');
disp('Figure 2 shows the results for the reference network');
disp('Please <Press a key> to continue estimation at all points');
if (do_profile~=1) pause; end;




%     根据最邻近参考网点来估计所有其他点
% --- Estimate all other points relatively to the closest ref. network point ---
disp(' ');
disp('Estimate the DEM error and displ. rate at all points.');
allESTPS    = zeros(NPM, NPS);% container for all estimated parameters
allvarf     = zeros(1,NPS);% all variance factors
invQy1      = inv(Qy1);% speed-up, outside loop for variance factor.
TRI         = DelaunayTri(X_REF',Y_REF');% easily find closest point
for p=1:NPS
  % --- Report progress ---------------------------------
  if (mod(p-1,100)==0)
    disp(['  ILS: point ', num2str(p), ' of ', num2str(NPS)]);
  end
  % --- Find the closest reference network point --------
  this_X      = X(p);
  this_Y      = Y(p);
  qrypts      =[this_X, this_Y];
  IDX_closest = nearestNeighbor(TRI,qrypts);
  X_closest   = X_REF(IDX_closest);
  Y_closest   = Y_REF(IDX_closest);
  IDX_this_refnetpnt = IDX_REF(IDX_closest);% "global" index
  ref_topo    = PM_at_REF(IDX_closest,1);%
  ref_defo    = PM_at_REF(IDX_closest,2);%
  % --- Check if it is a point that needs to be connected --------------
  dist_to_net = sqrt((this_X-X_closest).^2+(this_Y-Y_closest).^2);
  if (IDX_this_refnetpnt~=p)
    % --- 0: Double-difference phase observations at arc ---------------
    % --- dPhi=to-from --> dH=H(to)-H(from) --> H(to)=H(from)+dH
    y1 = wrap(PHASE(:,p) - PHASE(:,IDX_this_refnetpnt));
    % --- 1: float solution for ambiguities ----------------------------
    afloat = y1./(-2.*pi);% float solution (special case for PSI)
    zfloat = Z.' * afloat;% decorrelate float solution (for faster search)
    % --- 2: Obtain bound for ILS search, function "ebs" ---------------
    [zfixed_bs, sqnorm_bs]   = ebs(zfloat,L,D,1,invLtDL);% 1 candidate
    % --- 3: Fixed solution for decorrelated ambiguities using "ils"----
    [zfixed_ils, sqnorm_ils, ierr] = ils(zfloat,Linv,Dinv,sqnorm_bs,1);
    % --- Check ierr status flag ---------------------------------------
    if (bitand(ierr,1))
      warning('not enough candiates found in "ils"');
    end
    if (bitand(ierr,2))
      %warning('maxloop reached in "ils", using bootstrapped estimate');
      sqnorm_ils = sqnorm_bs;
      zfixed_ils = zfixed_bs;
    end
    % --- 4: Fixed solution a_check, Inverse Z-transform ---------------
    afixed = invZt*zfixed_ils;
    % --- 5: Unwrap data using the estimated ambiguities ---------------
    y_uw   = y1 + 2.*pi.*afixed;
    % --- 6: Ordinary weighted least-squares using unwrapped data ------
    bhat   = PROJ_LS*y_uw;% estimated DEM error, defo differences
    yhat   = B1*bhat;
    ehat   = wrap(y_uw-yhat);% residuals are always wrapped, closest solution
    allvarf(p)    = (ehat.'*invQy1*ehat)./(NIFG-NPM);
    allESTPS(:,p) = [ref_topo; ref_defo] + bhat;
  % --------------------------------------------------------------------
  % --- This is a point in the reference network -----------------------
  else
    if (dist_to_net>eps)
      warning('unexpected');
    end
    % --- At reference points, set the parameters directly -------------
    allvarf(p)    = -1;%
    allESTPS(:,p) = [ref_topo; ref_defo];
  end
end
disp(' ...done');
disp(' ');


%      实践中无法获取的真值参考点
% --- "true" at reference point not known in practice ------------------
topo_refpnt = simtopo(IDX_REF(IDX_refpnt));% "global" index
defo_refpnt = simdefo(IDX_REF(IDX_refpnt));% "global" index
% --- Report the errors ------------------------------------------------
errtopo     = simtopo-topo_refpnt - allESTPS(1,:);
errdefo     = simdefo-defo_refpnt - allESTPS(2,:);
disp('Comparing all estimates with simulated input:');
disp(['  Max. error DEM error:   ', num2str(max(abs(errtopo))), ' [m]']);
disp(['  Max. error Displ. rate: ', num2str(max(abs(errdefo))), ' [mm/y]']);

%     好的参考阈值
% --- Threshold good estimates -----------------------------------------
disp(' ');
disp('Thresholding on a posteriori variance factor');
IDX_ok = find(allvarf<2.0);
disp(['  Number of selected points: ', ...
  num2str(length(IDX_ok)), ' of ', num2str(length(errtopo))]);
disp(['  Mean error DEM error accepted points [m]:      ', ...
  num2str(mean(errtopo(IDX_ok)))]);
disp(['  Std. error DEM error accepted points [m]:      ', ...
  num2str(std(errtopo(IDX_ok)))]);
disp(['  Max. error DEM error accepted points [m]:      ', ...
  num2str(max(abs(errtopo(IDX_ok))))]);
disp(['  Mean error Displ. rate accepted points [mm/y]: ', ...
  num2str(mean(errdefo(IDX_ok)))]);
disp(['  Std. error Displ. rate accepted points [mm/y]: ', ...
  num2str(std(errdefo(IDX_ok)))]);
disp(['  Max. error Displ. rate accepted points [mm/y]: ', ...
  num2str(max(abs(errdefo(IDX_ok))))]);


%     绘制每个点的估计参数
% --- Plot the estimated parameters at the points ----------------------
figure(3)
  subplot(3,2,1)
    plotps(X, Y, simtopo-topo_refpnt, [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Simulated DEM errors')
  subplot(3,2,2)
    plotps(X, Y, simdefo-defo_refpnt, [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Simulated Displ. rates')
  subplot(3,2,3)
    plotps(X(IDX_ok), Y(IDX_ok), allESTPS(1,IDX_ok), [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Estimated DEM error [m] at PS points')
  subplot(3,2,4)
    plotps(X(IDX_ok), Y(IDX_ok), allESTPS(2,IDX_ok), [-40,40]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Estimated Displ. rates [mm/y] at PS points')
  subplot(3,2,5)
    plotps(X(IDX_ok), Y(IDX_ok), errtopo(IDX_ok), [-5,5]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Error with simulated input [m]')
  subplot(3,2,6)
    plotps(X(IDX_ok), Y(IDX_ok), errdefo(IDX_ok), [-5,5]);
    hold on
    plot(X_REF(IDX_refpnt), Y_REF(IDX_refpnt), 'r+', 'MarkerSize',8);
    hold off
    title('Error with simulated input [mm/y]')



% --- Save results ----------------------------------------------
disp('------------------------------------');
disp(' ');
key = [];
while (isempty(key))
  key = input('Enter "s" to save results, "q" to quit: ', 's');
end
if (key=='s') 
  [FILENAME, PATHNAME] = uiputfile('stundemo_data.mat', 'Save as');
  file = [PATHNAME, FILENAME];
  save(file, 'NIFG', 'NPS', 'maxDEM', 'stdDEFO', 'stdNSE', 'maxAPS', ...
    'Fk_true', 'X', 'Y', 'Bperp', 'Btemp', 'acq_times', ...
    'PHASE', 'PHASE_UW', 'simtopo', 'simdefo', 'phi_noise', 'phi_aps');
  disp(['Variables saved in ', FILENAME]);
  disp('You can load these variables to repeat processing.');
end



% --- End -------------------------------------------------------
disp(' ');
disp('End of demonstration');
disp(' ');
if (do_profile==1)
  figure(4);
  %profile report htmlfile
  profile plot;
  profile off;
end
more on
%%%EOF
