%-------------------------------------------------------------------------%
% CODEX_calc.m
% CODEX matrix calculation for a regular polygon of n spins
% Noah H. Somberg, Westley W. Wu, Mei Hong
% Written in MATLAB R2021b
% July 2022
% This script takes a distance matrix (which here is generated based on a
% regular polygon) and calculates a time-dependent CODEX decay using
% 1H-driven spin diffusion theory, and plots the results\
% See publication: https://doi.org/10.1021/acs.biochem.2c00464
%-------------------------------------------------------------------------%


%--------SIMULATION PARAMETERS--------%
gamma = 251.185e6; % Gyromagnetic ratio of fluorine
mu_0 = 1.25663706212e-6; % Vacuum permeability
hbar = 1.054571817e-34; % Planck constant
gammaProt = 267.52218744e6; % Proton gyro ratio
ang = 1e-10; % One angstrom
powd = 0.2; % Powder average of angular dependence

uplimit = 5000; % Upper bound of CODEX plot in ms (x axis)
step = 1; % Time increment in ms (smaller equals smoother curves)

n = 5; % Oligamer number
%s = 8; % Nearest neighbor distance

time_ax = 0:step:uplimit;
F0 = 3.4; % Overlap integra
s = 8.8; % NN distance in A

% Create a polygon with n sides of length n at [0,0] with initial rotation
% 0 degrees
poly = createPoly(n,s,[0;0],0); 
dismatrix = zeros(n,n); % Initialize a matrix for all distances

    
for p1 = 1:n
    for p2 = 1:n % For each set of coordinates
        xdist = abs(poly(1,p2) - poly(1,p1)); % Get x dist
        ydist = abs(poly(2,p2) - poly(2,p1)); % Get y dist
        dismatrix(p1,p2) = sqrt(xdist^2 + ydist^2); % Get total dist
    end
end
    
% Calculate known couplings to double check parameter values are correct
prot_1a = (mu_0 * hbar * gammaProt^2)/(4*pi*ang^3); % 1A 1H dipolar coup
prot_1a_hz = prot_1a/(2*pi); % Convert to Hz, should be 120120 Hz
    
f_1a = (mu_0 * hbar * gamma^2)/(4*pi*ang^3); % 1 A F-F dip coup (in rads!)
    
dipcoup = f_1a;
        
M0matrix = eye(n); % Initial state is identity matrix

np=uplimit/step+1; % Number of points

W=dismatrix.^(-3)*dipcoup; % Homonuclear dipolar coupling strength
Wsqu=W.^2; % Coupling squared

for i=1:n % Detailed balance
    Wsqu(i,i)=0; % Zero the diagonal of the coupling matrix
    Wsqu_sums = sum(Wsqu,1); % Sum each column
    Wsqu(i,i) = -Wsqu_sums(i); % Diag set to neg sum
end    
            
% Calculate the exchange matrix K
K=0.5*pi*Wsqu*powd*F0/1000000;
% Mt is a 3d matrix: 
% first axis is the ending spin, 
% second axis is starting spin, 
% 3rd axis is time
Mt = zeros(n,n,np); 
prop = expm(step/1000*K); % Calculate the propegator
    
for currSpin = 1:n
    % For each spin,  
    % calc the dip coup matrix W and the dip coup square
    
    M0 = M0matrix(:,currSpin); % Extract the vector for init mag on spin
    currMat = expm(0/1000*K); % Calculate the intial SD matrix
    
    t_idx = 1; % Initialize a time index

    for t = 0:step:uplimit
        % For each time step, calculate exchange process
        Mt(currSpin,:,t_idx) = currMat*M0; % Calculate Mt
        currMat = prop*currMat; % Increment exchange matrix

        t_idx = t_idx + 1; % Increment time index
    end
end
    
Mt_avg = zeros(2,np); % Initialize a matrix for avg magnetization
Mt_avg(1,:) = 0:step:uplimit; % First row is time, second row is M(t) avg

% Calculate average over all spins
for t_idx = 1:np
    % Calc avg mag at each time pt
    Mt_avg(2,t_idx) = trace(Mt(:,:,t_idx))/n; 
end

sim = Mt_avg(2,:);

% Plot the result
figure;
set(gca, 'FontName', 'Arial')
hold on
plot(time_ax./1000,sim,'k','LineWidth',2)
xlim([0,5]);
ylim([0, 1.1]);
xticks([0:1:5]);
yticks([0:0.2:1]);
box on
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',16)
set(gca,'linewidth',2)

function points = createPoly(n,s,origin,rot)
    % Outputs a set of points forming a regular polygon with n sides of
    % length s centered at origin, with an initial rotation with respect to
    % the x axis of rot

    theta = 360/n; % Angle for drawing polygon radial vectors
    points = zeros(2,n); % Initialize output matrix
    Rotmat = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    % Length of first vector based on side, using law of cosines
    veclen = sqrt(s^2/(2-2*cosd(theta))); 
    
    xtran = origin(1);
    ytran = origin(2);

    % Get first vertex from length of radial vector and initial incination
    % angle
    points(:,1) = [veclen*cosd(rot); veclen*sind(rot)]; 

    % Calculate remaining points by applying rotation matrix
    for i = 2:n
        points(:,i) = Rotmat*points(:,i-1); 
    end

    % Translate the shape to the specified origin
    points(1,:) = points(1,:)+xtran;
    points(2,:) = points(2,:)+ytran;
end

%-------------------------------------------------------------------------%
% End of CODEX_calc.m
%-------------------------------------------------------------------------%