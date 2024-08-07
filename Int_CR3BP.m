function [dX] = Int_CR3BP(t,X,prms)
% For numerical integration in the normalized CR3BP
% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u)

% Preallocate state output
dX = zeros(6,1);

% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+prms.u)^2 + y^2 + z^2);
r2 = sqrt((x+prms.u-1)^2 + y^2 + z^2);

% Equations of Motion - CR3BP
ddx = 2*dy + x - (1-prms.u)*(x+prms.u)/(r1^3) - prms.u*(x+prms.u-1)/(r2^3);
ddy = -2*dx + y -((1-prms.u)/(r1^3) + prms.u/(r2^3))*y;
ddz = -((1-prms.u)/(r1^3) + prms.u/(r2^3))*z;

% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end