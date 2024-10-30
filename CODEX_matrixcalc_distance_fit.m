%-------------------------------------------------------------------------%
% CODEX_matrixcalc_distance_fit.m
% CODEX matrix calculation for a regular polygon of n spins, while fitting
% the nearest-neighbor distance to the data points
% Noah H. Somberg, Westley W. Wu, Mei Hong
% Written in MATLAB R2021b in July 2022
% Revised October 2024
% See publication: https://doi.org/10.1021/acs.biochem.2c00464
%-------------------------------------------------------------------------%


%--------SIMULATION PARAMETERS--------%
close all
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
F0 = 3.41;


timePoints = [100 250 500 1000 1500 2000 3000 4000];
dataPoints = [0.74 0.56 0.49 0.35 0.29 0.24 0.19 0.14]; % In ms
err = [0.02 0.02 0.02 0.02 0.02 0.03 0.02 0.04];

r_ax = 4:0.1:15; % Axis to fit on (in Angstrom)
r_idx = 1;

% Initialize a matrix to store all sims
allCurves = zeros(length(time_ax),length(r_ax));


for s = r_ax
    disp(s)
    %--- Create oligomer distance matrix ---%
    olig = createPoly(n,s,[0;0],0); % Create one oligamer at 0,0    
    dismatrix = zeros(n,n); % Initialize a matrix for all distances
    
    allPoints = olig;
    
    for p1 = 1:n
        for p2 = 1:n % For each set of coordinates
            xdist = abs(allPoints(1,p2) - allPoints(1,p1)); % Get x dist
            ydist = abs(allPoints(2,p2) - allPoints(2,p1)); % Get y dist
            dismatrix(p1,p2) = sqrt(xdist^2 + ydist^2); % Get total dist
        end
    end
    
    
    %--- Initial Calculations ---%
    
    prot_1a = (mu_0 * hbar * gammaProt^2)/(4*pi*ang^3); % Calculate the proton 1 angstrom dipolar coupling, as a sanity check
    prot_1a_hz = prot_1a/(2*pi); % Convert to Hz, should be 120120 (i.e. 120.120 kHz)
    
    f_1a = (mu_0 * hbar * gamma^2)/(4*pi*ang^3); % Calculate 1 angstrom F-F dipolar coupling (in radians!)
    
    dipcoup = f_1a;
        
    M0matrix = eye(n); % Initial state is identity matrix
    
    np=uplimit/step+1; % Number of points
    
    W=dismatrix.^(-3)*dipcoup; % Homonuclear dipolar coupling strength
    Wsqu=W.^2; % Coupling squared
    
    for i=1:n % Detailed balance
        Wsqu(i,i)=0; % Zero the diagonal of the coupling matrix
        Wsqu_sums = sum(Wsqu,1); % Sum each column
        Wsqu(i,i) = -Wsqu_sums(i); % Set each diagonal to the negative sum of the column
    end    
            
    % Calculate the exchange matrix K
    K=0.5*pi*Wsqu*powd*F0/1000000;  % Ignore (1-3*costheta^2) by replacing it with 0.8.
    Mt = zeros(n,n,np); % Mt is a 3d matrix: first axis is the ending spin, second axis is starting spin, 3rd axis is time
    
    prop = expm(step/1000*K); % Calculate the propegator
    
    for currSpin = 1:n
        % For each spin, calculate the 
        %calculate the dipolar coupling matrix W and the dipolar coupling square .
        
        M0 = M0matrix(:,currSpin); % Extract the vector corresponding to the initial mag on this spin
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
    
    for t_idx = 1:np
        Mt_avg(2,t_idx) = trace(Mt(:,:,t_idx))/n; % Calc avg mag at each time point
    end
    
    allCurves(:,r_idx) = Mt_avg(2,:);
    r_idx = r_idx + 1;
end


chiSq = zeros(1,length(r_ax));
for r_idx = 1:length(r_ax)
    sumChi = 0;
    for timeIdx = 1:length(timePoints)
            % For each measured time point
            meas = dataPoints(timeIdx);
            predIdx = timePoints(timeIdx)+1;
            pred = allCurves(predIdx,r_idx); % Get the predicted value at the time point
            chiSqCalc =  (pred - meas)^2/err(timeIdx)^2;
            sumChi = sumChi + chiSqCalc;
            chiSq(r_idx) = sqrt(sumChi);
    end
end

[bestRMSD,bestIdx] = min(chiSq);

figure;
plot(r_ax,chiSq,'k','LineWidth',2);
yticks([0:10:50]);
ylim([0,50]);
xlim([4 15]);
box on
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',16)
set(gca,'linewidth',2)
xlabel('Nearest-Neighbor Distance (Anstrom)')
ylabel('Chi Squared Nu')

figure;
set(gca, 'FontName', 'Arial')
hold on
plot(time_ax./1000,allCurves(:,bestIdx),'k','LineWidth',2)
errorbar(timePoints./1000,dataPoints,err,'o','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',6,'MarkerFaceColor','k','CapSize',12,'Color','k')
xlim([0,5]);
ylim([0, 1.1]);
xticks([0:1:5]);
yticks([0:0.2:1]);
xlabel('Mixing Time (s)')
ylabel('S/S0')
box on
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',16)
set(gca,'linewidth',2)


disp('Best r')
disp(r_ax(bestIdx))
disp('Chi sq')
disp(chiSq(bestIdx))
disp('Red chi sq')
disp(chiSq(bestIdx)/(length(dataPoints)-1))