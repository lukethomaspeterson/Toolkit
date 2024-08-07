function [tcnf_input] = GenTCNFin(k,I,theta)

% Generate NF coords for TCNF input 
% Inputs:
%          k - Number of points in each action/angle pair [3x1]
%          I - List of action variables [3xN]
%          theta - List of intervals of angle variables [3x2]

I1 = I.I1; I2 = I.I2; I3 = I.I3; 

k1 = k(1); k2 = k(2); k3 = k(3);
N1 = length(I1); N2 = length(I2); N3 = length(I3); 

theta1 = linspace(theta(1,1),theta(1,2),k1);
theta2 = linspace(theta(2,1),theta(2,2),k2);
theta3 = linspace(theta(3,1),theta(3,2),k3);

temp = zeros(k1*k2*k3*N1*N2*N3,6); 
id = 1; 

for i1 = 1:N1
    for j1 = 1:k1
        for i2 = 1:N2
            for j2 = 1:k2
                for i3 = 1:N3
                    for j3 = 1:k3
                        temp(id,:) = [I1(i1), theta1(j1), I2(i2), theta2(j2), I3(i3), theta3(j3)]; 
                        id = id+1;
                    end
                end
            end
        end
    end
end

% Change coordinates 
tcnf_input = zeros(k1*k2*k3*N1*N2*N3,6); 
for i = 1:k1*k2*k3*N1*N2*N3
    tcnf_input(i,:) = [sqrt(temp(i,1)).*exp(temp(i,2)), sqrt(temp(i,1)).*exp(-temp(i,2)), sqrt(2*temp(i,3)).*cos(temp(i,4)), -sqrt(2*temp(i,3)).*sin(temp(i,4)), sqrt(2*temp(i,5)).*cos(temp(i,6)), -sqrt(2*temp(i,5)).*sin(temp(i,6))];
end

end