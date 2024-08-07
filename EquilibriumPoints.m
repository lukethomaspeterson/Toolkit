function r = EquilibriumPoints(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EquilibriumPoints computes the equilibrium points for the CR3BP given a %
% mass ratio between the two primary bodies mu = u                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  u      = mass ratio [1x1]                                       %
% Output: r      = (stacked) position vectors for the eq. pts. [5x3]      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(5,1);
y = x;
z = x; 

% Partials for Newton Method
dUdx = @(x) x - (1-u)*(x+u)/abs(x+u)^3 - u*(x-1+u)/abs(x-1+u)^3; % y=z=0 
d2Udx2 = @(x) 1 + 2*(1-u)/abs(x+u)^3 + 2*u/abs(x-1+u)^3; % y=z=0 

%%% Collinear equilibrium points: Use Newton Method 
% L1, Between P1 & P2 
x(1) = NewtonMethod(dUdx, d2Udx2, 1/2-u);
y(1) = 0;

% L2, Right of P2
x(2) = NewtonMethod(dUdx, d2Udx2, (1-u)+1e-5); 
y(2) = 0;

% L3, Left of P1
x(3) = NewtonMethod(dUdx, d2Udx2, 1.25*(-u)); 
y(3) = 0;

%%% Triangular equilibrium points: Know explicitly 
% L4, Upper Triangular Eq. Pt.
x(4) = 1/2 - u;
y(4) = sqrt(3)/2;

% L5, Lower Triangular Eq. Pt. 
x(5) = 1/2 - u;
y(5) = -sqrt(3)/2;

r = [x,y,z]; 

end