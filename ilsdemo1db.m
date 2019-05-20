%ILSDEMO1DB   Basic demonstration of ILS estimation
%   This script generates a 1d polynomial plus noise 
%     y = b1*x + b2*x^2 + e,
%   wraps it, and the ILS estimator is used to estimate
%   the parameters b1 and b2.
%   
%   These observations y are assumed to be uncorrelated,
%   i.e., a diagonal vc-matrix is used.
%   
%   Random data are generated each time this script is 
%   executed.
%
%   Example:
%     ilsdemo1db
%
%   See also STUN, ILSDEMO1D.

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
% $Revision: 1.7 $  $Date: 2005/12/21 10:31:08 $

more off
disp('ILSDEMO1DB   Basic demonstration of ILS estimation');
disp('   This script generates a 1d polynomial plus noise');
disp('     y = b1*x + b2*x^2 + e,');
disp('   wraps it, and the ILS estimator is used to estimate');
disp('   the parameters b1 and b2.');
disp(' ');

while (1)
% --- Remove the reference results if they exist ---------
if (exist('my_sqnorm_bs')) clear my_*; end


% --- Load or simulate the data --------------------------
disp(' ');
disp('===========================================================');
disp('Enter "1" to load a good example to check your installation');
disp('      "2" to load an example gone bad');
disp('      "q" to quit');
disp('      "-" to generate random input');
disp('===========================================================');
disp(' ');
key = input('Please enter: ','s');
if (isempty(key)) key = '-'; end;
switch (lower(key))
  case '1',
    if ~exist('poly_good.mat')
      error('cannot find file "poly_good.mat"');
    end
    load poly_good
    disp('Loading data set "poly_good.mat"');
    disp('This is an example to demonstrate a successful estimation.');
    disp('The reference results are stored in the variables:');
    whos my_*
    % After run, the results were save in poly_good.mat by:
    % y = y1;
    % my_sqnorm_bs  = sqnorm_bs;
    % my_sqnorm_ils = sqnorm_ils;
    % my_zfixed_ils = zfixed_ils;
    % my_zfixed_bs  = zfixed_bs;
    % my_afloat     = afloat;
    % my_afixed_ils = afixed_ils;
    % save poly_good N b1 b2 b x nl n y a my_*
  case '2',
    if ~exist('poly_wrong.mat')
      error('cannot find file "poly_wrong.mat"');
    end
    load poly_wrong
    disp('Loading data set "poly_wrong.mat"');
    disp('For this set, there exists an incorrect solution that has a');
    disp('large 2nd degree term, such that the unwrapped model fits');
    disp('better than the correct model.');
    disp('If the variance of the second pseudo-observation');
    disp('could be constrained more, this may have been avoided.');
    disp(' ');
    %previously variables were saved:
    %save poly_wrong N b1 b2 b x nl n y a
  case 'q',
    more on;
    return;
  otherwise,
    % --- Simulate the line and the noise --------------------
    rand('state',sum(100*clock));% change state of random number generator
    randn('state',sum(100*clock));% change state of random number generator
    N  = 5+round(25*rand);% number of samples in [5,30]
    b1 = 5.*randn;% true slope of line;
    b2 = 2.*randn;% bit smaller 2nd order term
    b  = [b1;b2];% "b" contains the true float parameters
    x  = 10.*rand(N,1)-5;% times of samples in [-5,5]
    nl = 20.*pi./180 + 1.0*rand;% std of noise on data [rad]
    n  = nl.*randn(N,1);% realization of noise
    [y, a] = wrap([x,x.^2]*b+n);% "y" simulated observations; "a" true ambiguities
end

disp(['Number of samples: ', num2str(N)]);
disp(['Noise level [deg]: ', num2str(nl.*180./pi)]);
disp(['Slope of line:     ', num2str(b1)]);
disp(['Second coefficient:', num2str(b2)]);
disp(['Number of wraps:   ', num2str(max(a)-min(a))]);
disp(' ');
disp('Press a key to view the simulated data');
disp(' ');
pause



% --- Plot the simulated data ----------------------------
xx  = linspace(-5,5,400).';% continuous axis for plots
fig = figure(1);
set(fig,'Position',[100 200 800 400]);
h = plot(xx,[xx,xx.^2]*b,'b-', ...
     xx,wrap([xx,xx.^2]*b),'b--', ...
     x,y,'r.');
set(h(1),'LineWidth',3);
title('Wrapped observations: y = b1*x + b2*x^2 - 2*pi*a');
xlabel('x (time in years)')
ylabel('y (phase in radians)')



% --- Set up model and estimate ILS -----------------------
disp(' ');
disp('Setting up the model of observation equations:');
disp('   y1 = A1*a + B1*b + e1;  Qy1=D{e1},              (Eq. 1)');
disp('where:');
disp('  y1 [N x 1] vector with N wrapped phase observations.');
disp('  A1 [N x N] design matrix for integer parameters.');
disp('  a  [N x 1] vector of N unknown integer parameters (the ambiguities).');
disp('  B1 [N x 1] design matrix for float parameters.');
disp('  b  [1 x 1] vector of 2 unknown float parameters.');
disp('  e1 [N x 1] vector of measurement noise.');
disp(' Qy1 [N x N] vc-matrix of the noise (D{e} is dispersion of e).');
disp(' ');
disp('This system of equations (Eq. 1) is under-determined,');
disp('i.e., there are more unknown parameters than observations.');
disp('Therefor, (Eq. 1) is regularized using 2 pseudo-observation:');
disp('  [y1] = [A1] * a + [B1] * b + [e1];  Qy=[D{e1}  0 ],  (Eq. 2)');
disp('  [y2]   [A2]       [B2]       [e2]      [ 0  D{e2}]');
disp(' ');

% --- Add pseudo-observation -------------------------------
NPM = length(b);% number of float parameters == number of pseudo-obs.
y1  = y;
y2  = zeros(NPM,1);% zero pseudo-observations
y   = [y1; y2];

% --- Set up design matrix for integer parameters ----------
A1  = -2*pi*eye(N);
A2  = zeros(NPM,N);
A   = [A1; A2];

% --- Set up design matrix for float parameters ------------
B1  = [x, x.^2];% functional model
B2  = eye(NPM);
B   = [B1; B2];

% --- Set up vc-matrix -------------------------------------
var_obs        = 1^2;% assumed a priori variance on observations
var_pseudo_obs = [10^2, 5^2];% "soft-bounds" on parameters
Qy1 = var_obs*eye(N);
Qy2 = diag(var_pseudo_obs);
Qy  = [Qy1,     zeros(N,NPM); ...
       zeros(NPM,N),    Qy2];

% --- Wait for user ----------------------------------------
disp('Set up of matrices completed');
disp(' ');
disp('Press a key to estimate the parameters using ILS');
disp(' ');
pause



% --- Perform the estimation using the LAMBDA method -------
% --- 1: Float solution for integer parameters -------------
% --- Alternatively the float solution for the integer -----
% --- parameters can be computed using the solution for a --
% --- partioned system, Eq 3.15 ----------------------------
P_B    = eye(N+NPM) - B * inv(B.'*inv(Qy)*B) * B.'*inv(Qy);
A_     = P_B*A;% reduced design matrix
Qahat  = inv(A_.'*inv(Qy)*A_);%
%afloat = Qahat*A_.'*inv(Qy)*y;% float solution partioned model
afloat = y1./(-2.*pi);% float solution (special case for PS)


% --- 2: Decorrelation using Z-transform, function "zt" ----
[Z,L,D,ZtQyZ,zfloat]         = zt(Qahat,afloat);


% --- 3: Obtain bound for ILS search, function "ebs" -------
disp('Searching for best two fits');
ncands = 2;% get the two best fitting slopes
[zfixed_bs, sqnorm_bs]   = ebs(zfloat,L,D,ncands);


% --- 4: Fixed solution decorrelated ambiguities using "ils"-
Chi2 = max(sqnorm_bs);
Linv = inv(L);
Dinv = 1./D;
[zfixed_ils, sqnorm_ils, ierr] = ils(zfloat,Linv,Dinv,Chi2,ncands);
% --- Check ierr status flag -------------------------------
if (bitand(ierr,1))
  warning('not enough candiates found in "ils"');
end
if (bitand(ierr,2))
  warning('maxloop reached in "ils", using bootstrapped estimate');
  zfixed_ils = zfixed_bs;
end


% --- 5: Fixed solution acheck, Inverse Z-transform --------
afixed_ils = round(inv(Z.'))*zfixed_ils;% best two solutions
disp('Estimates obtained using ILS (and bootstrap)');
disp(' ');


% --- 6: Unwrap data using the estimated ambiguities -------
y_uw     = y1 + 2.*pi.*afixed_ils(:,1);%best solution
y_uw_2nd = y1 + 2.*pi.*afixed_ils(:,2);%second best solution


% --- 7: Estimate float parameter using unwrapped data -----
% --- Least-squares solution of E{y}=Ax; D{y}=Qy -----------
% --- is given by xhat=inv(A.'*inv(Qy)*A)*A.'*inv(Qy)*y ----
% --- We use this direct experession here for clarity ------
% --- However, there are faster and more stable ways -------
bhat     = inv(B1.'*inv(Qy1)*B1)*B1.'*inv(Qy1) * y_uw;
bhat_2nd = inv(B1.'*inv(Qy1)*B1)*B1.'*inv(Qy1) * y_uw_2nd;

% --- Compute the adjusted residuals and the a -------------
% --- posteriori variance factor ---------------------------
yhat       = B1*bhat;% adjusted observations
ehat       = wrap(y_uw-yhat);% adjusted residuals
redundancy = N-NPM;% 2 parameters estimated using N observations
varf       = (ehat.'*inv(Qy1)*ehat) ./ redundancy;
nl_post    = sqrt(varf*var_obs)*180./pi;%noise level

% ---
yhat_2nd   = B1*bhat_2nd;% adjusted observations
ehat_2nd   = wrap(y_uw-yhat_2nd);% adjusted residuals
coh_2nd    = enscoh(ehat_2nd);


% --- Report results ---------------------------------------
disp('Estimated slope, best fit');
disp('--------------------------------');
disp(['  number of samples:                  ', num2str(N)]);
disp(['  noise level [deg]:                  ', num2str(nl.*180./pi)]);
disp(['  simulated parameters b:             ', num2str(b.')]);
disp(['  estimated parameters bhat:          ', num2str(bhat.')]);
disp(['  error with input:                   ', num2str((b-bhat).')]);
disp(['  a priori variance factor:           ', num2str(var_obs)]);
disp(['  a posteriori variance factor:       ', num2str(varf)]);
disp(['  --> estimated noise level [deg]:    ', num2str(nl_post)]);
disp(['  ensemble coherence:                 ', num2str(enscoh(ehat))]);
err_amb = round(sum(abs(a-afixed_ils(:,1))));
disp(['  summed error estimated ambiguities: ', num2str(err_amb)]);
if (err_amb==0)
  disp( '  ambiguities estimated correctly:    yes     <----');
else
  disp( '  ambiguities estimated correctly:    no      <---- **');
end
disp(' ');
disp('Estimated slope, second best fit');
disp('--------------------------------');
disp(['  second best fitting slope bhat:     ', num2str(bhat_2nd.')]);
disp(['  ensemble coherence:                 ', num2str(enscoh(ehat_2nd))]);
err_amb = round(sum(abs(a-afixed_ils(:,2))));
disp(['  summed error estimated ambiguities: ', num2str(err_amb)]);
disp(' ');
disp('Press a key to plot results');
disp(' ');
pause



% --- Plot results on top of simulated data -----------------
figure(fig)
hold on;
BBB = [xx,xx.^2];% functional model
plot(xx,BBB*bhat,'g-', ...  % estimated best slope
     x,y_uw,'gx');
plot(xx,BBB*bhat_2nd,'y-', ...  % estimated second best slope
     x,y_uw_2nd,'y+');
hold off;
legend('unwrapped model', ...
       'wrapped model', ...
       'observed wrapped phase', ...
       'estimated model', ...
       'estimated unwrapped phase', ...
       '2nd best estimate', -1);
       %'2nd best estimate unwrapped phase');

disp('End of demonstration');
disp('Run it again for new random input simulation');
disp(' ');


% --- Check with my results if they were loaded from disk --------------
if (exist('my_sqnorm_bs'))
  if (max(abs(my_sqnorm_bs-sqnorm_bs))>1d-3 | ...
      max(abs(my_sqnorm_ils-sqnorm_ils))>1d-3 | ...
      max(abs(my_zfixed_bs-zfixed_bs))>1d-3 | ...
      max(abs(my_zfixed_ils-zfixed_ils))>1d-3 | ...
      max(abs(my_afloat-afloat))>1d-3 | ...
      max(abs(my_afixed_ils-afixed_ils))>1d-3)
    error('difference detected between your results and the reference results!')
  else
    disp('Your results are identical to the reference results!');
  end
end
end %while loop
more on
%%%EOF
