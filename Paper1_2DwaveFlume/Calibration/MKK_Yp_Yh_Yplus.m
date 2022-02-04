%% Calculate Yp from idealized or realized y+ values
%CFD modeling reques an estimate of Yp, the height of smallest grid cell
%centroid near the wall boundary. The smallest grid cell height is
%therefore 2*Yp for a given y+ value. This script also estimates the
%maximum cell size from the desired number of grid refinement levels, and
%compares the non-dimensional bed roughness coefficient (Ks+) with y+ to
%determine if using a roughness approach is appropriate for the model. 
%
% This is version 1 of this script.
%
% BKN - USGS PCMSC 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate simulation Reynolds Number:
rho = 1000; %density of fluid used in model
U = 0.05; %characteristic velocity (use depth-averaged from model!)
L = 1.9; %characteristic length scale (e.g., water depth, coral height?)
nu = 1E-6; %kinematic viscosity used in model
yPlus = 16; %either a "guess" or calculated from the model
Ks = 123E-6*30; %Nikuradse roughness length (estimated from velocity profiles)
numRefinementLevels = 5; %number of refinement levels for max cell height

mu = nu*rho;

Re = (rho*U*L)/mu;

%Calculate skin friction coefficient
cf = (2*log10(Re)-6.5)^-2.3;

%Calculate wall shear stress
tauw = (0.5*rho*(U^2))*cf;

%Calculate friction velocity at the wall
ut = sqrt(tauw/rho);

%Calculate Yp
Yp = (yPlus*mu)/(rho*ut);
fprintf('Yp = %0.4f m\n',Yp);

%Calculate Yh (height of first cell)
fprintf('min Cell height = %0.4f m\n',2*Yp);

%Calculate maximum cell height
fprintf('max Cell height = %0.4f m\n',2*Yp*(2^numRefinementLevels));

%Calculate Ks+ (dimensionless roughness height)
KsPlus = (rho*Ks*ut)/mu;
fprintf('Ks+ = %0.2f\n',KsPlus)
fprintf('y+ = %0.2f\n',yPlus)
fprintf('Note: Ks+ should be smaller than y+ for roughness to be contained\nwithin the first grid cell near the wall\n')
fprintf('If Ks+ >> 90, using a roughness approach at the wall is not appropriate!\n')