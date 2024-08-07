function [QP] = ReadOut(XdX)

% Generate NF coords from Cartesian state
% Input: XdX   - Cartesian Position/Velocity [kx6]
% Output: QP   - Position/Momenta pairs for NF change of variables [kx6]

% CoV: Cartesian to Hamiltonian
% XdX = (x,y,z,dx,dy,dz) 
% CoV: Step 1 -- (x,y,z,dx,dy,dz) --> (qx,px,qy,py,qz,pz)
QP_US = zeros(size(XdX));
QP_US(:,[1 3 5 6]) = XdX(:,[1 2 3 6]); % position & dz coordinates
QP_US(:,2) = XdX(:,4) - XdX(:,2); % px = dx - y
QP_US(:,4) = XdX(:,5) + XdX(:,1); % py = dy + x

% CoV: Step 2 -- USian to Spanish model transformation (+ --> -)
QP = -QP_US; 

end