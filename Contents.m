% STUN algorithm for Persistent Scatterers Interferometry Toolbox.
%
% Elementary Functions.
%   bs_success  - Success rate of the simple bootstrap estimator.
%   enscoh      - Ensemble coherence.
%   plotarc     - Plot arcs of network in color.
%   plotps      - Plot PS points in color.
%   sparsify    - Bin points in grid cells.
%   wrap        - Wrap phase data.
%
% Functional Model.
%   simacq      - Simulate acquisition baselines.
%   simphi      - Simulate ERS-like phase observations.
%   simpos      - Simulate 2D positions.
%
% Stochastic Model.
%   psivce      - Variance Component Estimation for PSI.
%   psivcmtx    - VC-matrix for double-difference observations.
%
% Ambiguity Resolution.
%   ebs         - Extended bootstrap fixed solution.
%   ils         - Integer least-squares fixed solution.
%   ltdl        - LTDL decomposition Q=L.'*D*L.
%   zt          - Z-transformation (decorrelation).
%
% Demonstrations.
%   ilsdemo1d   - ILS estimation of the slope of wrapped line.
%   ilsdemo1db  - ILS estimation of 2nd degree polynomial.
%   stundemo    - Main demonstration of STUN algorithm.
%   vcedemo     - Variance Component Estimation for PSI.
%
% Data Sets.
%   poly_good.mat    - Example data set for ilsdemo1db.
%   poly_wrong.mat   - Example data set for ilsdemo1db.
%   vcedemodat.mat   - Example data set/reference results (vcedemo).
%   stundemodat1.mat - Fractal displacement example data (stundemo).
%   stundemodat2.mat - Fractal atmosphere example data (stundemo).
%


% The STUN toolbox accompanies the book:
% Kampes, B.M., "Radar Interferometry -- The Persistent
% Scatterer Technique", published by Springer, 2006.

% You are allowed to use and/or modify this code for your
% own purposes.

% The STUN toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Written by (c)Bert Kampes, 21-Oct-2005
% Tested with Matlab version: 5.3.0.10183 (R11)
% $Revision: 1.6 $  $Date: 2005/12/21 10:31:08 $

