function H = Hamil(I,nf) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hamil.m reads in a list of action variables I = (I1,I2,I3) & the table  %
% of coefficients from nf.res, then returns the value of the truncated    %
% Hamiltonian function evaluated at the point(s) I = (I1,I2,I3).          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:                                                                 %
% .. I  - Nx3 vector of actions (I1,I2,I3)                                %
% .. nf - Mx5 matrix with exponents in cols 1-3, coefficients in col 4    %
%         col 5 is imaginary parts of coefficients & are ignored          %
%                                                                         %
% Output:                                                                 %
% .. H - N-vector of the evaluated truncated Hamiltonian                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Separate actions 
I1 = I(:,1); % saddle, commonly all 0s 
I2 = I(:,2); % center 
I3 = I(:,3); % center 

% Separate nf.res file into exponents & coefficients 
k1 = nf(:,1); % I1^k1 
k2 = nf(:,2); % I2^k2 
k3 = nf(:,3); % I3^k3
c = nf(:,4); % coefficients 

N = length(I1); % number of points (I1,I2,I3) 
M = length(c); % number of coefficients 

% Pre-allocate 
temp = 0; 
H = zeros(N,1);

for i = 1:length(I1)
    for j = 1:M
        temp = temp + c(j)*I1(i)^k1(j)*I2(i)^k2(j)*I3(i)^k3(j);
    end
    H(i) = temp; 
    temp = 0; 
end 
