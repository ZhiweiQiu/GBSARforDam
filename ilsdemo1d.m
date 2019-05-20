%ILSDEMO1D   Basic demonstration of ILS estimation
%   This script generates a line plus noise 
%     y = b*x + e,
%   wraps it, and the Integer Least-Squares (ILS) 
%   estimator is used to estimate the slope b.
%   
%   The observations y are assumed to be uncorrelated,
%   i.e., a diagonal vc-matrix is used.
%   
%   Random data are generated each time this script is 
%   executed.
%
%   Example:
%     ilsdemo1d
%
%   See also STUN, ILSDEMO1DB.

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

more off
disp('ILSDEMO1D   Basic demonstration of ILS estimation');
disp('   This script generates a line plus noise y = b*x + e,');
disp('   wraps it, and the ILS estimator is used to estimate');
disp('   the slope b.');
disp(' ');
disp('Press a key to generate random input');
disp(' ');
pause


% --- Simulate the line and the noise --------------------
rand('state',sum(100*clock));% change state of random number generator
randn('state',sum(100*clock));% change state of random number generator
N  = 5+round(25*rand);% number of samples in [5,30]
b  = 5.*randn;% slope of line; "b" is the true float parameter
x  = 10.*rand(N,1)-5;% times of samples in [-5,5]
nl = 20.*pi./180 + 1.0*rand;% std of noise on data [rad]
n  = nl.*randn(N,1);% realization of noise
[y, a] = wrap(b*x+n);% "y" simulated observations; "a" true ambiguities
disp(['Number of samples: ', num2str(N)]);
disp(['Noise level [deg]: ', num2str(nl.*180./pi)]);
disp(['Slope of line:     ', num2str(b)]);
disp(['Number of wraps:   ', num2str(max(a)-min(a))]);
disp(' ');
disp('Press a key to view the simulated data');
disp(' ');
pause



% --- Plot the simulated data ----------------------------
xx  = linspace(-5,5,400).';% continuous axis for plots
fig = figure(1);
set(fig,'Position',[100 200 800 400]);
plot(xx,b*xx,'b-', ...
     xx,wrap(b*xx),'b--', ...
     x,y,'r.');
title('Wrapped observations: y = x*b - 2*pi*a');
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
disp('  b  [1 x 1] vector of 1 unknown float parameters (the slope).');
disp('  e1 [N x 1] vector of measurement noise.');
disp(' Qy1 [N x N] vc-matrix of the noise (D{e} is dispersion of e).');
disp(' ');
disp('This system of equations (Eq. 1) is under-determined,');
disp('i.e., there are more unknown parameters than observations.');
disp('Therefor, (Eq. 1) is regularized using 1 pseudo-observation,');
disp('for the slope, with value 0 and standard deviation 10.0');
disp('  [y1] = [A1] * a + [B1] * b + [e1];  Qy=[D{e1}  0 ],  (Eq. 2)');
disp('  [y2]   [A2]       [B2]       [e2]      [ 0  D{e2}]');
disp(' ');

% --- Add pseudo-observation -------------------------------
NPM = 1;% number of float parameters == number of pseudo-obs.
y1  = y;
y2  = zeros(NPM,1);% zero pseudo-observation(s)
y   = [y1; y2];

% --- Set up design matrix for integer parameters ----------
A1  = -2*pi*eye(N);
A2  = zeros(NPM,N);
A   = [A1; A2];

% --- Set up design matrix for float parameters ------------
B1  = x;% x*slope
B2  = eye(NPM);
B   = [B1; B2];

% --- Set up vc-matrix -------------------------------------
var_obs        = 1^2;% assumed a priori variance on observations
var_pseudo_obs = 10^2;% assumed "soft-bound" on slope
Qy1 = var_obs*eye(N);
Qy2 = var_pseudo_obs*eye(NPM);
Qy  = [Qy1,     zeros(N,NPM); ...
       zeros(NPM,N),    Qy2];

% --- Wait for user ----------------------------------------
disp('Set up of matrices completed');
disp(' ');
disp('Press a key to estimate the slope using ILS');
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
[zfixed_bs, zsqnorm_bs]   = ebs(zfloat,L,D,ncands);


% --- 4: Fixed solution decorrelated ambiguities using "ils"-
Chi2 = max(zsqnorm_bs);
Linv = inv(L);
Dinv = 1./D;
[zfixed_ils, zsqnorm_ils, ierr] = ils(zfloat,Linv,Dinv,Chi2,ncands);
% --- Check ierr status flag -------------------------------
if (bitand(ierr,1))
  warning('not enough candiates found in "ils"');
end
if (bitand(ierr,2))
  warning('maxloop reached in "ils", using bootstrapped estimate');
  zfixed_ils = zfixed_bs;
end


% --- 5: Fixed solution acheck, Inverse Z-transform --------
afixed_ils = inv(Z.')*zfixed_ils;% best two solutions
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
yhat       = bhat*x;% adjusted observations
ehat       = wrap(y_uw-yhat);% adjusted residuals
redundancy = N-1;% 1 parameter estimated using N observations
varf       = (ehat.'*inv(Qy1)*ehat) ./ redundancy;



% --- Report results ---------------------------------------
disp('Estimated slope, best fit');
disp('--------------------------------');
disp(['  number of samples:                  ', num2str(N)]);
disp(['  noise level [deg]:                  ', num2str(nl.*180./pi)]);
disp(['  simulated slope:                    ', num2str(b)]);
disp(['  estimated slope bhat:               ', num2str(bhat)]);
disp(['  error with input:                   ', num2str(b-bhat)]);
disp(['  a priori variance factor:           ', num2str(var_obs)]);
disp(['  a posteriori variance factor:       ', num2str(varf)]);
nl_post = sqrt(varf*var_obs)*180./pi;
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
disp(['  second best fitting slope bhat:     ', num2str(bhat_2nd)]);
coh = enscoh(wrap(y_uw_2nd-bhat_2nd*x));
disp(['  ensemble coherence:                 ', num2str(coh)]);
err_amb = round(sum(abs(a-afixed_ils(:,2))));
disp(['  summed error estimated ambiguities: ', num2str(err_amb)]);
disp(' ');
disp('Press a key to plot results');
disp(' ');
pause



% --- Plot results on top of simulated data -----------------
figure(fig)
hold on;
plot(xx,xx*bhat,'g-', ...  % estimated best slope
     x,y_uw,'gx');
plot(xx,xx*bhat_2nd,'y-', ...  % estimated second best slope
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
more on
%%%EOF
