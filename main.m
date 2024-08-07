%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR3BP Local Orbital Elements Toolkit - main script                      %
%                                                                         %
% In this script, the user will be able to:                               %
% - Choose a Lagrange points of interest, Earth-Moon L1 or L2             %
% - Choose the order of normal form series expansion, if desired...       %
%   - If no specified order chosen, automatically sets to 10th order      %
% - Transform coordinates between Cartesian and local orbital elements    %
% - Generate periodic and quasi-periodic orbits from local action-angle   %
%   orbital elements -- note: orbits in the normal form dynamics          %
% - Compute stable and unstable manifolds of orbits selected via          %
%   action-angle elements                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; % clear workspace


%% User-input data

% NOTE: If selecting L2, go into Toolkit\data and replace all cvnf.1-6 with
%       Toolkit\data\L2\cvnf.1-6 by copying from the folder. 

% Select Lagrange Point of interest
LagrangePoint = input("Select Earth-Moon Lagrange Point of interest, L1 or L2. [Enter '1' for L1, '2' for L2]\n");

% If no Lagrange Point is selected, default to L1
if isempty(LagrangePoint)
    LagrangePoint = 1;
end

% Load normal form data based on user selection of Lagrange Point
if LagrangePoint == 1
    nf = dlmread('data\L1\nf.res');
elseif LagrangePoint == 2
    nf = dlmread('data\L2\nf.res');
end

ExpansionOrderActions = sum(nf(end,1:3));
ExpansionOrder = 2*ExpansionOrderActions;

% System data: Earth-Moon CR3BP 

% Universal Gravitational Constant 
G = 6.6743015e-20; % km3 kg-1 s-1

% Earth
bodies.earth.name  = 'earth';
bodies.earth.u     = 398600.435436; % km^3 * s^-2
bodies.earth.color = [235,214,173]./255;
bodies.earth.mass  = bodies.earth.u/G; % kg
bodies.earth.a     = 149598023; % km
bodies.earth.e     = 0.016708617; % [-]
bodies.earth.R     = 6371; % km, volumetric mean radius 

% Moon 
bodies.moon.name    = 'moon';
bodies.moon.u       = 4902.800066; % km^3 / s^2
bodies.moon.color   = [0, 1, 1];
bodies.moon.mass    = bodies.moon.u/G; % kg
bodies.moon.a       = 384400; % km
bodies.moon.e       = 0.05490; % [-] 
bodies.moon.R       = 1737.4; % km, volumetric mean radius 
bodies.moon.R_n     = bodies.moon.R / bodies.moon.a; % normalized radius
bodies.moon.MR      = bodies.moon.mass / (bodies.moon.mass + bodies.earth.mass); % Mass ratio w/ primary
bodies.moon.MR      = 0.012150584269542001;

% Compute Earth normalized radius in Earth-Moon frame
bodies.earth.R_n    = bodies.earth.R / bodies.moon.a; % normalized radius 

% Choose bodies
primary   = bodies.earth; 
secondary = bodies.moon; 

% Normalizing constants
lNorm = secondary.a; % n <-> km
tNorm = sqrt(lNorm^3/(primary.u+secondary.u)); % n <-> sec
vNorm = lNorm / tNorm; % n <-> km/sec

% Shortcut variables
mu   = secondary.MR;
R1_n = primary.R_n;
R2_n = secondary.R_n;

% Equilibrium Points
rLPs = EquilibriumPoints(mu);

% Lagrange Points
L1 = rLPs(1,:);
L2 = rLPs(2,:);

% Plot Primary and Secondary Bodies
[Xsphere,Ysphere,Zsphere] = sphere; % Generate sphere
sx1 = Xsphere*primary.R_n - mu; % Earth
sy1 = Ysphere*primary.R_n; % Earth
sz1 = Zsphere*primary.R_n; % Earth
sx2 = Xsphere*R2_n + (1-mu); % Moon
sy2 = Ysphere*R2_n; % Moon
sz2 = Zsphere*R2_n; % Moon

% Colors for Primary and Secondary Bodies
colors.ltblue = [125 216 255]./255;
colors.grey   = [127 127 127]./255;

% ------------------------------------------------- 
% Integration Options
% -------------------------------------------------
tol               = 1e-13; % integration toloerance
options           = odeset('RelTol',tol,'AbsTol',tol); % set ODE integration options
options_Primary   = odeset('Event',@event_Primary,'RelTol',tol,'AbsTol',tol);
options_Secondary = odeset('Event',@event_Secondary,'RelTol',tol,'AbsTol',tol);
prms.u            = mu; % save mass ratio into parameter structure
prms.R2_n         = R2_n; % save lunar radius into parameter structure


%% Example 1: Generate Energy Level Sets

% Example configured for:
% - Lagrange Point: L1
% - Expansion Order: 8

% Set number of points in grid (I = (I1,I2,I3))
N = 5e2;

% Generate points
if LagrangePoint == 1
    if (ExpansionOrderActions < 16) && (ExpansionOrderActions > 10)
        I2_max = 0.4;
        I3_max = 0.5;
    elseif ExpansionOrderActions == 10
        I2_max = 0.5;
        I3_max = 0.6;
    elseif ExpansionOrderActions < 10
        I2_max = 0.6;
        I3_max = 0.7;
    end
elseif LagrangePoint == 2
    if (ExpansionOrderActions < 16) && (ExpansionOrderActions > 10)
        I2_max = 0.6;
        I3_max = 0.7;
    elseif ExpansionOrderActions == 10
        I2_max = 0.7;
        I3_max = 0.8;
    elseif ExpansionOrderActions < 10
        I2_max = 0.9;
        I3_max = 1;
    end
end
[I2,I3] = meshgrid(linspace(0,I2_max,N),linspace(0,I3_max,N)); % for contours
I = zeros(N^2,3); % initialize (I2,I3) space
I(:,2) = reshape(I2,N^2,1); % set up I2 for contour plot
I(:,3) = reshape(I3,N^2,1); % set up I3 for contour plot

% Compute H
H = Hamil(I,nf);

% Prepare for contour plot
Hplot = reshape(H,N,N); % reshape H

% Plot contours
v0 = (2:2:20)./10; % initialize contour values
v  = repelem(v0,2); % modify for contour function (needs repeated values)

figure
[C,h] = contour(I2,I3,Hplot,v);
clabel(C,h,v)
axis square
xlabel('$I_2$','interpreter','latex')
ylabel('$I_3$','interpreter','latex')
if LagrangePoint == 1
    title('Energy Level Sets of $H$ about EM $L_1$','Interpreter','latex')
elseif LagrangePoint == 2
    title('Energy Level Sets of $H$ about EM $L_2$','Interpreter','latex')
end


%% Example 2a: Generate Periodic Orbit

clear I

% Choose orbit amplitude, i.e., planar action variable I2
I.I1 = 0;
I.I2 = 0.1;
I.I3 = 0;

% Generate angles parameterizing orbit
n_PO = 500; % number of points along PO
k = [1;n_PO;1]; % k - Number of points in each action/angle pair [3x1]
theta = [0 0;0 2*pi;0 2*pi]; % List of intervals of angle variables [3x2]

% Transform to NF coordinates
tcnf_input = GenTCNFin(k,I,theta);

% Change variables to Cartesian state
% ..Set up .txt file
if LagrangePoint == 1
    fileID = fopen('L1\PlanarLyapunov.txt','w');
elseif LagrangePoint == 2
    fileID = fopen('L2\PlanarLyapunov.txt','w');
end
stri = [num2str(length(tcnf_input(:,1))) ' ' num2str(length(tcnf_input(1,:)))];
fprintf(fileID,stri)
fclose(fileID);
if LagrangePoint == 1
    writematrix(tcnf_input,'L1\PlanarLyapunov.txt','Delimiter',' ','WriteMode','append')
elseif LagrangePoint == 2
    writematrix(tcnf_input,'L2\PlanarLyapunov.txt','Delimiter',' ','WriteMode','append')
end

% ..Compute change of variables
if LagrangePoint == 1
! .\bin\tcnf.exe ..\L1\PlanarLyapunov.txt ..\L1\PlanarLyapunov_syn.txt 1
elseif LagrangePoint == 2
! .\bin\tcnf.exe ..\L2\PlanarLyapunov.txt ..\L2\PlanarLyapunov_syn.txt 1
end

% Read-in synodic Cartesian coordinates
if LagrangePoint == 1
    PO_syn = ReadIn('L1\PlanarLyapunov_syn.txt');
elseif LagrangePoint == 2
    PO_syn = ReadIn('L2\PlanarLyapunov_syn.txt');
end

% Plot orbit in synodic rotating frame
figure();
hold on
plot(rLPs(LagrangePoint,1),rLPs(LagrangePoint,2),'rd','MarkerFaceColor','r')
plot(PO_syn(:,1),PO_syn(:,2),'k','LineWidth',1)
axis equal
grid on
xlabel('$x$ [-]','Interpreter','latex')
ylabel('$y$ [-]','Interpreter','latex')
title('Planar Lyapunov Orbit','Interpreter','latex')
if LagrangePoint == 1
    legend('$L_1$',strcat('$I_2$ = ',num2str(I.I2(1))),'interpreter','latex')
elseif LagrangePoint == 2
    legend('$L_2$',strcat('$I_2$ = ',num2str(I.I2(1))),'interpreter','latex')
end


%% Example 3: Generate Stable/Unstable Manifolds for a Periodic Orbit
% Function(s): GenTCNFin.m, WsuNFseeds.m, ReadIn.m 

clear I

% Choose orbit amplitude, i.e., planar action variable I2
I.I1 = 0;
I.I2 = 0.1;
I.I3 = 0;

% Generate angles parameterizing orbit
n_PO = 500; % number of points along PO
k = [1;n_PO;1]; % k - Number of points in each action/angle pair [3x1]
theta = [0 0;0 2*pi;0 2*pi]; % List of intervals of angle variables [3x2]

% Transform to NF coordinates
tcnf_input = GenTCNFin(k,I,theta);

% Generate data for coordinate transformation to NF coordinates
Wsu_input = WsuNFseeds(tcnf_input);

% Change variables to Cartesian state
% ..Set up .txt file
if LagrangePoint == 1
    fileID = fopen('L1\PlanarLyapunov_Wsu.txt','w');
elseif LagrangePoint == 2
    fileID = fopen('L2\PlanarLyapunov_Wsu.txt','w');
end
stri = [num2str(length(Wsu_input(:,1))) ' ' num2str(length(Wsu_input(1,:)))];
fprintf(fileID,stri)
fclose(fileID);
if LagrangePoint == 1
    writematrix(Wsu_input,'L1\PlanarLyapunov_Wsu.txt','Delimiter',' ','WriteMode','append')
elseif LagrangePoint == 2
    writematrix(Wsu_input,'L2\PlanarLyapunov_Wsu.txt','Delimiter',' ','WriteMode','append')
end

% ..Compute change of variables
if LagrangePoint == 1
! .\bin\tcnf.exe ..\L1\PlanarLyapunov_Wsu.txt ..\L1\PlanarLyapunov_Wsu_syn.txt 1
elseif LagrangePoint == 2
! .\bin\tcnf.exe ..\L2\PlanarLyapunov_Wsu.txt ..\L2\PlanarLyapunov_Wsu_syn.txt 1
end

% Read-in synodic Cartesian coordinates
if LagrangePoint == 1
    PO_Wsu_syn = ReadIn('L1\PlanarLyapunov_Wsu_syn.txt');
elseif LagrangePoint == 2
    PO_Wsu_syn = ReadIn('L2\PlanarLyapunov_Wsu_syn.txt');
end

% Save each Wsu seed separately 
xUp0 = PO_Wsu_syn(1:n_PO,:);
xUn0 = PO_Wsu_syn((n_PO+1):(2*n_PO),:);
xSp0 = PO_Wsu_syn((2*n_PO+1):(3*n_PO),:);
xSn0 = PO_Wsu_syn((3*n_PO+1):(4*n_PO),:);

% Integration time
t0 = 0;
tf = 10;

% Compute the stable/unstable manifolds of the orbit 
figure
hold on
plot(rLPs(LagrangePoint,1),rLPs(LagrangePoint,2),'rd','MarkerFaceColor','r')
surf(sx2,sy2,sz2,'FaceColor',colors.grey,'EdgeColor',colors.grey); view(2)
surf(sx1,sy1,sz1,'FaceColor',colors.ltblue,'EdgeColor',colors.ltblue); view(2)
if LagrangePoint == 1
    for i = 1:5:length(xSp0)
        % Integrate stable/unstable manifolds
        [tSp, XSp] = ode113(@Int_CR3BP, [tf t0], xSp0(i,:), options_Primary, prms);
        [tSn, XSn] = ode113(@Int_CR3BP, [tf t0], xSn0(i,:), options_Secondary, prms);
        [tUp, XUp] = ode113(@Int_CR3BP, [t0 tf], xUp0(i,:), options_Secondary, prms);
        [tUn, XUn] = ode113(@Int_CR3BP, [t0 tf], xUn0(i,:), options_Primary, prms);
    
        % Plot the stable/unstable manifolds in blue/red, respectively
        plot(XSp(:,1),XSp(:,2),'b','LineWidth',1) % Stable "positive" --> Primary
        plot(XUp(:,1),XUp(:,2),'r','LineWidth',1) % Unstable "positive" --> Secondary
        plot(XSn(:,1),XSn(:,2),'b','LineWidth',1) % Stable "negative" --> Secondary
        plot(XUn(:,1),XUn(:,2),'r','LineWidth',1) % Unstable "negative" --> Primary
    end
elseif LagrangePoint == 2
    for i = 1:5:length(xSp0)
        % Integrate stable/unstable manifolds
        [tSp, XSp] = ode113(@Int_CR3BP, [tf t0], xSp0(i,:), options_Secondary, prms);
        [tSn, XSn] = ode113(@Int_CR3BP, [tf t0], xSn0(i,:), options_Primary, prms);
        [tUp, XUp] = ode113(@Int_CR3BP, [t0 tf], xUp0(i,:), options_Primary, prms);
        [tUn, XUn] = ode113(@Int_CR3BP, [t0 tf], xUn0(i,:), options_Secondary, prms);
    
        % Plot the stable/unstable manifolds in blue/red, respectively
        plot(XSp(:,1),XSp(:,2),'b','LineWidth',1) % Stable "positive" --> Primary
        plot(XUp(:,1),XUp(:,2),'r','LineWidth',1) % Unstable "positive" --> Secondary
        plot(XSn(:,1),XSn(:,2),'b','LineWidth',1) % Stable "negative" --> Secondary
        plot(XUn(:,1),XUn(:,2),'r','LineWidth',1) % Unstable "negative" --> Primary
    end
end
plot(PO_syn(:,1),PO_syn(:,2),'k','LineWidth',1)
if LagrangePoint == 1
legend('$L_1$','Moon','Earth','Stable','Unstable','interpreter','latex')
elseif LagrangePoint == 2
legend('$L_2$','Moon','Earth','Stable','Unstable','interpreter','latex')
end
xlabel('$x$ [-]','interpreter','latex')
ylabel('$y$ [-]','interpreter','latex')
if LagrangePoint == 2
xlim([0.91 1.69])
ylim([-0.3 0.3])
end
axis equal
grid on 


%% Example 4: Generate Periodic Orbits

clear I

% Choose orbit amplitude, i.e., planar action variable I2
I.I1 = 0;
I.I2 = 0.05:0.05:0.6;
I.I3 = 0;

% Generate angles parameterizing orbit
n_PO = 200; % number of points along PO
k = [1;n_PO;1]; % k - Number of points in each action/angle pair [3x1]
theta = [0 0;0 2*pi;0 2*pi]; % List of intervals of angle variables [3x2]

% Transform to NF coordinates
tcnf_input = GenTCNFin(k,I,theta);

% Change variables to Cartesian state
% ..Set up .txt file
if LagrangePoint == 1
    fileID = fopen('L1\PlanarLyapunovFamily.txt','w');
elseif LagrangePoint == 2
    fileID = fopen('L2\PlanarLyapunovFamily.txt','w');
end
stri = [num2str(length(tcnf_input(:,1))) ' ' num2str(length(tcnf_input(1,:)))];
fprintf(fileID,stri)
fclose(fileID);
if LagrangePoint == 1
    writematrix(tcnf_input,'L1\PlanarLyapunovFamily.txt','Delimiter',' ','WriteMode','append')
elseif LagrangePoint == 2
    writematrix(tcnf_input,'L2\PlanarLyapunovFamily.txt','Delimiter',' ','WriteMode','append')
end

% ..Compute change of variables
if LagrangePoint == 1
! .\bin\tcnf.exe ..\L1\PlanarLyapunovFamily.txt ..\L1\PlanarLyapunovFamily_syn.txt 1
elseif LagrangePoint == 2
! .\bin\tcnf.exe ..\L2\PlanarLyapunovFamily.txt ..\L2\PlanarLyapunovFamily_syn.txt 1
end

% Read-in synodic Cartesian coordinates
if LagrangePoint == 1
    PO_fam_syn = ReadIn('L1\PlanarLyapunovFamily_syn.txt');
elseif LagrangePoint == 2
    PO_fam_syn = ReadIn('L2\PlanarLyapunovFamily_syn.txt');
end

c = parula(length(I.I2));

figure
hold on 
plot(rLPs(LagrangePoint,1),rLPs(LagrangePoint,2),'rd','MarkerFaceColor','r')
surf(sx2,sy2,sz2,'FaceColor',colors.grey,'EdgeColor',colors.grey); view(2)
for i = 1:length(I.I2)
    plot(PO_fam_syn((i-1)*n_PO+1:i*n_PO,1),PO_fam_syn((i-1)*n_PO+1:i*n_PO,2),'Color',c(i,:))
end
hold off
axis equal 
grid on 
xlabel('$x$ [-]','Interpreter','latex')
ylabel('$y$ [-]','Interpreter','latex')
if LagrangePoint == 1
legend('$L_1$','Moon','interpreter','latex')
elseif LagrangePoint == 2
legend('$L_2$','Moon','interpreter','latex')
end
colorbar;
clim([min(I.I2), max(I.I2)+0.0001]);
title('Planar Lyapunov Orbits','interpreter','latex')


%% Example 5: Generate Quasi-Periodic Orbits
% Function(s): GenTCNFin.m

clear I c

% Choose orbit amplitude, i.e., planar action variable I2
I.I1 = 0;
I.I2 = 0.025:0.05:0.175; % define planar actions / amplitudes
I.I3 = 0.1; % define vertical actions / amplitudes

% Generate angles parameterizing orbit
n_QPO = 50; % number of points along PO
N_QPO = length(I.I2)*length(I.I3); % number of QPOs computed
k = [1;n_QPO;n_QPO]; % k - Number of points in each action/angle pair [3x1]
theta = [0 0;0 2*pi;0 2*pi]; % List of intervals of angle variables [3x2]

% Transform to NF coordinates
tcnf_input = GenTCNFin(k,I,theta);

% Change variables to Cartesian state
% ..Set up .txt file
if LagrangePoint == 1
    fileID = fopen('L1\LissajousFamily.txt','w');
elseif LagrangePoint == 2
    fileID = fopen('L2\LissajousFamily.txt','w');
end
stri = [num2str(length(tcnf_input(:,1))) ' ' num2str(length(tcnf_input(1,:)))];
fprintf(fileID,stri)
fclose(fileID);
if LagrangePoint == 1
    writematrix(tcnf_input,'L1\LissajousFamily.txt','Delimiter',' ','WriteMode','append')
elseif LagrangePoint == 2
    writematrix(tcnf_input,'L2\LissajousFamily.txt','Delimiter',' ','WriteMode','append')
end

% ..Compute change of variables
if LagrangePoint == 1
! .\bin\tcnf.exe ..\L1\LissajousFamily.txt ..\L1\LissajousFamily_syn.txt 1
elseif LagrangePoint == 2
! .\bin\tcnf.exe ..\L2\LissajousFamily.txt ..\L2\LissajousFamily_syn.txt 1
end

% Read-in synodic Cartesian coordinates
if LagrangePoint == 1
    QPO_fam_syn = ReadIn('L1\LissajousFamily_syn.txt');
elseif LagrangePoint == 2
    QPO_fam_syn = ReadIn('L2\LissajousFamily_syn.txt');
end

% Save each orbit
K = k(1)*k(2)*k(3);
QPOs = zeros(K,6,1);
for i = 1:N_QPO
    QPOs(:,:,i) = QPO_fam_syn((i-1)*K+1:i*K,:);
end

% Plot Lissajous orbit as a surface
Xs = zeros(n_QPO,n_QPO+1,N_QPO);
Ys = zeros(n_QPO,n_QPO+1,N_QPO);
Zs = zeros(n_QPO,n_QPO+1,N_QPO);
for i = 1:N_QPO
    Xs(1:n_QPO,1:n_QPO,i) = reshape(QPOs(:,1,i),n_QPO,n_QPO,1)';
    Xs(:,end,i) = Xs(1:n_QPO,1,i);
    Ys(1:n_QPO,1:n_QPO,i) = reshape(QPOs(:,2,i),n_QPO,n_QPO,1)';
    Ys(:,end,i) = Ys(1:n_QPO,1,i);
    Zs(1:n_QPO,1:n_QPO,i) = reshape(QPOs(:,3,i),n_QPO,n_QPO,1)';
    Zs(:,end,i) = Zs(1:n_QPO,1,i);
end

c = copper(N_QPO);

figure
plot3(rLPs(LagrangePoint,1),rLPs(LagrangePoint,2),rLPs(LagrangePoint,3),'rd','MarkerFaceColor','r')
hold on
for i = 1:N_QPO
    surf(Xs(:,:,i),Ys(:,:,i),Zs(:,:,i),'FaceColor',c(i,:),'FaceAlpha',0.5,'EdgeColor','none');
end
hold off
axis equal
grid on
rotate3d on
xlabel('$x$ [-]','interpreter','latex') 
ylabel('$y$ [-]','interpreter','latex')
zlabel('$z$ [-]','interpreter','latex')
title('Lissajous Orbits: $I_3 = 0.1$','interpreter','latex')
colormap(copper)
cb = colorbar;
clim([min(I.I2), max(I.I2)+0.0001]);
view((-1)^(LagrangePoint+1)*45,25) % flip viewing angle depending on Lagrange Point


%% Example 6: Local Orbital Elements of Numerically-Corrected Periodic Orbit

% Using pre-computed periodic orbit data files
% - Orbits considered: L1 & L2 Planar Lyapunov Orbits

% Load orbit data file
if LagrangePoint == 1
    % load L1 orbit
    OrbitData = load('L1\NumPO_syn.mat');
    PeriodicOrbit = OrbitData.XPO_L1(:,1:6);
    t             = OrbitData.XPO_L1(:,end);
    T             = t(end);
elseif LagrangePoint == 2
    % load L2 orbit
    OrbitData = load('L2\NumPO_syn.mat');
    PeriodicOrbit = OrbitData.XPO_L2(:,1:6);
    t             = OrbitData.XPO_L2(:,end);
    T             = t(end);
end

% Prepare for normal form coordinate transformation
% - change to position / momenta
tcnf_input = ReadOut(PeriodicOrbit);

% Change variables to local action-angle orbital elements
% ..Set up .txt file
if LagrangePoint == 1
    fileID = fopen('L1\NumPO_syn.txt','w');
elseif LagrangePoint == 2
    fileID = fopen('L2\NumPO_syn.txt','w');
end
stri = [num2str(length(tcnf_input(:,1))) ' ' num2str(length(tcnf_input(1,:)))];
fprintf(fileID,stri)
fclose(fileID);
if LagrangePoint == 1
    writematrix(tcnf_input,'L1\NumPO_syn.txt','Delimiter',' ','WriteMode','append')
elseif LagrangePoint == 2
    writematrix(tcnf_input,'L2\NumPO_syn.txt','Delimiter',' ','WriteMode','append')
end

% ..Compute change of variables
if LagrangePoint == 1
! .\bin\tcnf.exe ..\L1\NumPO_syn.txt ..\L1\NumPO.txt -1
elseif LagrangePoint == 2
! .\bin\tcnf.exe ..\L2\NumPO_syn.txt ..\L2\NumPO.txt -1
end

% Read-in local orbital elements
if LagrangePoint == 1
tcnf_output = dlmread('L1\NumPO.txt');
elseif LagrangePoint == 2
tcnf_output = dlmread('L2\NumPO.txt');
end

QP_NF = tcnf_output(2:end,1:6); % [q1 p1 q2 p2 q3 p3]

% Compute local action-angle elements
A = NaN(size(QP_NF));
A(:,1) = QP_NF(:,1).*QP_NF(:,2); % I1
A(:,2) = -log(sqrt(QP_NF(:,1))./sqrt(QP_NF(:,2))); % theta1
A(:,3) = 0.5*(QP_NF(:,3).^2 + QP_NF(:,4).^2); % I2
A(:,4) = -atan2(QP_NF(:,4),QP_NF(:,3)); % theta2
A(:,5) = 0.5*(QP_NF(:,5).^2 + QP_NF(:,6).^2); % I3
A(:,6) = -atan2(QP_NF(:,6),QP_NF(:,5)); % theta3

% Plot action-angle elements
FigurePosition = [250,250,1190,390];

f=figure();
f.Position = FigurePosition;
subplot(2,4,[1 2 5 6])
hold on
plot(rLPs(LagrangePoint,1),rLPs(LagrangePoint,2),'rd','MarkerFaceColor','r')
plot(PeriodicOrbit(:,1),PeriodicOrbit(:,2),'k','LineWidth',1)
hold off
axis equal
grid on
xlabel('$x$ [-]','interpreter','latex')
ylabel('$y$ [-]','interpreter','latex')
title('Numerically-corrected Periodic Orbit','interpreter','latex')

subplot(2,4,3)
plot(t,QP_NF(:,1),'k')
% semilogy(t,A(:,1),'k') % plot saddle action I1 instead
grid on
xlim([t(1),t(end)])
xlabel('$t$ [-]','interpreter','latex')
ylabel('$q_1$ [-]','interpreter','latex')
% ylabel('$I_1$ [-]','interpreter','latex') % label I1

subplot(2,4,7)
plot(t,QP_NF(:,2),'k')
% plot(t,A(:,2),'k') % plot saddle angle theta1 instead
grid on
xlim([t(1),t(end)])
xlabel('$t$ [-]','interpreter','latex')
ylabel('$p_1$ [-]','interpreter','latex')
% ylabel('$\theta_1$ [-]','interpreter','latex') % label theta1

subplot(2,4,4)
plot(t,A(:,3),'k')
grid on
xlim([t(1),t(end)])
xlabel('$t$ [-]','interpreter','latex')
ylabel('$I_2$ [-]','interpreter','latex')

subplot(2,4,8)
plot(t,wrapTo2Pi(A(:,4)),'k')
grid on
xlim([t(1),t(end)])
ylim([0 2*pi])
yticks([0 pi/2, pi, 3*pi/2 2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('$t$ [-]','interpreter','latex')
ylabel('$\theta_2$ [-]','interpreter','latex')

