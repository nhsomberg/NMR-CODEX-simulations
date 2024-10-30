%-------------------------------------------------------------------------%
% CODEX_genetic_asymmetric_fit.m
% Implements a very simple genetic algorithm to find a set of distances for
% the specified N-mer that satisfies the CODEX data.
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
step = 10; % Time increment in ms (smaller equals smoother curves)

n = 5; % Oligamer number
%s = 8; % Nearest neighbor distance

time_ax = 0:step:uplimit;
F0 = 3.41;

% Experimental data for fitting
timePoints = [100 250 500 1000 1500 2000 3000 4000];
dataPoints = [0.74 0.56 0.49 0.35 0.29 0.24 0.19 0.14];
err = [0.02 0.02 0.02 0.02 0.02 0.03 0.02 0.04];
S0 = [1.00 0.96 0.87 0.83 0.74 0.68 0.55 0.41];

% Nearest-neighbor distance (within one oligamer)
initialNN = 10; % Initial guess for nearest-neighbor distance
randMove = 0.1; % Space to take random walks in one generation
genSize = 100; % Size of each generation
genMax = 50;

genSeed = createPoly(n,initialNN,[0;0],0);

timeAx = 0:step:uplimit;
nTimePoints = length(timeAx);

bestChis = zeros(1,genMax);
bestSims = zeros(genMax,nTimePoints);

savedSeeds = zeros(2,n,genMax);

% Uncomment for ensemble plotting
% for i = 1:100
%     hold on
%     scatter(gen1(1,:,i),gen1(2,:,i));
% end
% hold off

for genNum = 1:genMax
    genSims = zeros(genSize,nTimePoints); % Initialize matrix to store sims for this generation

    savedSeeds(:,:,genNum) = genSeed;

    % Generate initial ensemble
    angles = rand(1,n,genSize) .* 360;
    deltaX = randMove.*cosd(angles);
    deltaY = randMove.*sind(angles);
    delta = [deltaX; deltaY];
    delta(:,:,1) = zeros(2,n); % First page of deltas is zeros (intial guess always tested)
    gen = repmat(genSeed,1,1,genSize); % Populate gen1 with copies of initial points;
    gen = gen + delta;

    for currSim = 1:genSize

        currPoints = gen(:,:,currSim); % Get the coordinate for current sim

        dismatrix = zeros(n,n); % Initialize distance matrix

        % Calculate distance matrix
        for p1 = 1:n
            for p2 = 1:n % For each set of coordinates
                xdist = abs(currPoints(1,p2) - currPoints(1,p1)); % Get x dist
                ydist = abs(currPoints(2,p2) - currPoints(2,p1)); % Get y dist
                dismatrix(p1,p2) = sqrt(xdist^2 + ydist^2); % Get total dist
            end
        end

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

        genSims(currSim,:) = Mt_avg(2,:); % Save sims for this generation
    end % End loop for one generation

    % Now fit each curve in the current generation
    genChis = zeros(genSize,1);
    for currSim = 1:genSize
        chisq = 0;

        sim = genSims(currSim,:); % Get the current simulation

        for timeIdx = 1:length(timePoints)
            % For each measured time point
            simTimeIdx = timePoints(timeIdx)/step + 1;
            pred = sim(simTimeIdx); % Get the predicted value at the time point
            chisq = chisq + (((dataPoints(timeIdx)-pred)^2)/(err(timeIdx)^2)); % Calculate chi squared
        end
        genChis(currSim) = chisq;
    end

    [minVal, minIdx] = min(genChis); % Get the best fit
    bestChis(genNum) = minVal; % Save the best chi
    bestSims(genNum,:) = genSims(minIdx,:);

    genSeed = gen(:,:,minIdx);
end


% for i = 1:genMax
%     hold on;
%     plot(timeAx,bestSims(i,:))
% end

% hold on
% plot(timeAx,bestSims(genMax,:),'r','LineWidth',1.5)
% errorbar(timePoints,ss0vals,err,'ko','MarkerFaceColor','k','MarkerSize',8,'CapSize',10,'LineWidth',1.5);
% ylim([0, 1]);
% xlim([0, 5000])
% xlabel("CODEX Mixing Time (ms)");
% ylabel("S/S_0")

figure;
set(gca, 'FontName', 'Arial')
hold on
plot(time_ax./1000,bestSims(genMax,:),'k','LineWidth',2)
%scatter(timePoints./1000,dataPoints,75,'ko','LineWidth',1.5,'MarkerFaceColor','w')
errorbar(timePoints./1000,dataPoints,err,'o','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',6,'MarkerFaceColor','k','CapSize',12,'Color','k')
%scatter(timePoints./1000,S0,36,'o','MarkerEdgeColor','k','LineWidth',2)
xlim([0,5]);
ylim([0, 1.1]);
xticks([0:1:5]);
yticks([0:0.2:1]);
box on
set(gca, 'FontName', 'Arial')
set(gca,'FontSize',16)
set(gca,'linewidth',2)
hold off

figure
plot(1:genMax,bestChis./2,'k','LineWidth',1.5);
H=gca;
H.FontName='Arial';
H.LineWidth=1; %change to the desired value
H.FontSize=16;
xlabel('Generation Number');
ylabel('Fit (\chi^2)');


figure
hold on
meanLen = mean(sqrt(sum((genSeed.^2)))); % Get mean radial length
avgBest = createPoly(n,meanLen,[0;0],0);
%scatter(genSeed(1,:),genSeed(2,:));
ck = ones(n,3);
patch(avgBest(1,:),avgBest(2,:),'k','FaceColor','none','LineWidth',2);
patch(genSeed(1,:),genSeed(2,:),'r','FaceColor','none','EdgeColor','r','LineWidth',2);
%set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
axis equal
xlim([-10 10])
ylim([-10 10])
box on


% for i = 1:n
%     finDist(i) = sqrt(abs(genSeed(1,i) - avgBest(1,i))^2 + abs(genSeed(2,i) - avgBest(2,i))^2)
% end