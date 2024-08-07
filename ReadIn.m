function [XdX] = ReadIn(filename)

% Generate NF coords for TCNF input
% Input:  filename - name of TCNF output file, .txt
% Output: XdX      - Cartesian state [kx6]

% Read in NF code text file 
tcnf_output = dlmread(filename); % 'tcnf_output.txt' 

% CoV: Hamiltonian to Cartesian 
qp_Sp = tcnf_output(2:end,1:6); % [q1 p1 q2 p2 q3 p3]

% CoV: Step 1 -- (q,p) --> (x,dx) 
xdx_Sp = zeros(size(qp_Sp));
xdx_Sp(:,[1 3 5 6]) = qp_Sp(:,[1 3 5 6]);
xdx_Sp(:,2) = qp_Sp(:,2) + qp_Sp(:,3); % dx = px + y (px = dx - y)
xdx_Sp(:,4) = qp_Sp(:,4) - qp_Sp(:,1); % dy = py - x (py = dy + x)

% CoV: Step 2 -- [x dx y dy z dz] --> [x y z dx dy dz]
XdX_Sp = [xdx_Sp(:,1) xdx_Sp(:,3) xdx_Sp(:,5) xdx_Sp(:,2) xdx_Sp(:,4) xdx_Sp(:,6)];

% CoV: Step 3 -- Spanish to USian model transformation 
%             -- 180* rotation about z-axis
XdX = zeros(size(XdX_Sp));
XdX(:,1) = -XdX_Sp(:,1); % x --> -x 
XdX(:,2) = -XdX_Sp(:,2); % y --> -y
XdX(:,3) = XdX_Sp(:,3); % z --> z
XdX(:,4) = -XdX_Sp(:,4); % dx --> -dx
XdX(:,5) = -XdX_Sp(:,5); % dy --> -dy
XdX(:,6) = XdX_Sp(:,6); % dz --> dz

end