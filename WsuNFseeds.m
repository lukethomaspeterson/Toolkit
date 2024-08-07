function Wsu = WsuNFseeds(BSS) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WsuNFseeds.m reads in a list of points corresponding to a bdd special   %
% solution, e.g. equilibrium point, periodic/quasi-periodic orbit, in     %
% normal form coordinates (Q,P). WsuNFseeds returns a list of points 4x   %
% as long, divided into 4 sections, q1 = +/- eps, then p1 = +/- eps.      %
% This function is used to generate initial points that, when transformed %
% into Cartesian coordinates, serve as initial X0 for integration &       %
% computation of stable/unstable manifolds.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:                                                                 %
% .. BSS - Nx6 list of points corresponding to a bounded special solution %
%                                                                         %
% Output:                                                                 %
% .. Wsu - 4Nx6 list of points corresponding to initial stable/unstable   %
%    manifold points for numerical integration (following transformation) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Size of perturbation 
eps = 1e-4;

% Get size of list 
N = length(BSS(:,1));

% Pre-allocate output 
Wsu = zeros(4*N,6);

for i = 1:N
    Wsu(i,:) = [eps,BSS(i,2:6)]; % q1=+eps, p1=0
    Wsu(N+i,:) = [-eps,BSS(i,2:6)]; % q1=-eps, p1=0
    Wsu(2*N+i,:) = [0,eps,BSS(i,3:6)]; % q1=0, p1=+eps
    Wsu(3*N+i,:) = [0,-eps,BSS(i,3:6)]; % q1=0, p1=-eps
end 