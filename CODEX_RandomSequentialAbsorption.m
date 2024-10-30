%-------------------------------------------------------------------------%
% CODEXRSA_RandomSequentialAbsorption.m
% Streamlined and optimized simulation of CODEX of multiple oligomers in a
% bilayer. Oligomers are added to the bilayer by Random Sequential
% Absorption (RSA)
% Written in MATLAB R2023b
% July 2024
% Written by Noah Somberg, contributions from Westley W. Wu, Mei Hong
% See publication: https://doi.org/10.1021/acs.biochem.2c00464
%-------------------------------------------------------------------------%

% USER PARAMETERS
w = 500; % Width of the simulated membrane patch in Angstrom
LP_ratio = 15; % Protein monomer to PHOSPHOLIPID ratio
lipid_area = 60; % Area per lipid molecule, in angstrom^2
olig_num = 5; % Oligomer number to simulate
pore_collision_radius = 10; % Radius of the pore, for which no other pore should be within!
CODEX_nn_dist = 7.4; % Nearest-neighbor FF-distance
max_sim_time = 4000; % Maximum simulation time in ms
step_sim_time = 10; % Simulation step time: smaller takes longer, but makes smoother curves
coupling_limit = 50; % (angstrom) The program ignores any F-F couplings that are further apart than this

num_pores = ceil(w^2 / (((LP_ratio/2)*olig_num*lipid_area)+(pi*pore_collision_radius^2))); % Calculate how many pores should be on the membrane patch
centers = RSA_circle_corrdinates(w,num_pores,pore_collision_radius,true); % Get the oligomer centers
spin_coords = centers_to_spins(olig_num,centers,CODEX_nn_dist,true); % Get the coordinates for all the spins

figure
scatter(spin_coords(:,1),spin_coords(:,2))


dis_matrix = pairwise_distances(spin_coords); % Calculate the distance matrix

output = CODEX(dis_matrix,max_sim_time,step_sim_time,coupling_limit); % CODEX calculation

figure
plot(0:step_sim_time:max_sim_time,output,'k-','LineWidth',1)
yline(1/olig_num)

% === FUNCTIONS BELOW HERE === 

function coords = RSA_circle_corrdinates(w,num_pores,pore_collision_radius,plotBool)
    
    pore_coordinates = [rand*w, rand*w];
    
    for pore_index = 2:num_pores
    
        x_candidate = rand*w;
        y_candidate = rand*w;
    
       % While any of the distances to the other pores are less than twice the
       % collisions radius, keep re-rolling coordinates
        while sum(sqrt(sum((pore_coordinates - [x_candidate,y_candidate]).^2,2)) < 2*pore_collision_radius)
            % While rejected, re-roll
            %disp("reject")
            x_candidate = rand*w;
            y_candidate = rand*w;
        end
        
        % Once accepted, add it in!
        pore_coordinates(pore_index,:) = [x_candidate,y_candidate];
    
    end

    coords = pore_coordinates;
    
    if plotBool
        figure(101)
        hold on
        for pore_index = 1:num_pores
            x_left = pore_coordinates(pore_index,1)-pore_collision_radius;
            y_bottom = pore_coordinates(pore_index,2)-pore_collision_radius;
            rectangle('Position',[x_left,y_bottom,2*pore_collision_radius,2*pore_collision_radius],'Curvature',[1 1])
        end
        axis square
    end
end

function points = createPoly(n,s,origin,rot)
    % Outputs a set of points forming a regular polygon with n sides of
    % length s centered at origin, with an initial rotation with respect to
    % the x axis of rot

    theta = 360/n; % Angle for drawing polygon radial vectors
    points = zeros(2,n); % Initialize output matrix
    Rotmat = [cosd(theta), -sind(theta); sind(theta), cosd(theta)]; % Rotation matrix
    veclen = sqrt(s^2/(2-2*cosd(theta))); % Length of first vector based on side, using law of cosines
    
    xtran = origin(1);
    ytran = origin(2);

    % Get first vertex from length of radial vector and initial incination
    % angle
    points(:,1) = [veclen*cosd(rot); veclen*sind(rot)]; 

    for i = 2:n
        points(:,i) = Rotmat*points(:,i-1); % Calculate remaining points by applying rotation matrix
    end

    % Translate the shape to the specified origin
    points(1,:) = points(1,:)+xtran;
    points(2,:) = points(2,:)+ytran;
end

function coords = centers_to_spins(olig_num,centers,CODEX_nn_dist,plotBool)
    num_pores = height(centers);
    coords = [];
    for index = 1:num_pores
        curr_center = centers(index,:);
        rotation = rand*360;
        olig_coords = createPoly(olig_num,CODEX_nn_dist,curr_center,rotation);

        coords = [coords; olig_coords'];

    end

    if plotBool
        figure(101)
        hold on
        for m_idx = 1:height(coords)
            m_x = coords(m_idx,1);
            m_y = coords(m_idx,2);

            ahel_rad = 2.3; % Alpha-helix radius

            x_left = m_x-ahel_rad;
            y_bottom = m_y-ahel_rad;
            rectangle('Position',[x_left,y_bottom,2*ahel_rad,2*ahel_rad],'Curvature',[1 1])
        end
    end
end

function distance_matrix = pairwise_distances(spin_coords)
    N = height(spin_coords);
    distance_matrix = zeros(N);

    for spin1 = 1:N
        for spin2 = 1:N
            distance_matrix(spin1,spin2) = norm(spin_coords(spin1,:) - spin_coords(spin2,:));
        end
    end

end

function curve_out = CODEX(distance_matrix, time_upperlimit, time_step, max_distance_computed)
    gamma = 251.185e6; % Gyromagnetic ratio of fluorine
    mu_0 = 1.25663706212e-6; % Vacuum permeability
    hbar = 1.054571817e-34; % Planck constant
    gammaProt = 267.52218744e6; % Proton gyro ratio
    ang = 1e-10; % One angstrom
    powd = 0.2; % Powder average of angular dependence
    F0 = 3.41; % Overlap integral
    
    time_axis = 0:time_step:time_upperlimit; % Time axis for CODEX
    N = height(distance_matrix);
    num_time_points = length(time_axis);

    % This is just a sanity check for units
    prot_1a = (mu_0 * hbar * gammaProt^2)/(4*pi*ang^3); % Calculate the proton 1 angstrom dipolar coupling, as a sanity check
    prot_1a_hz = prot_1a/(2*pi); % Convert to Hz, should be 120120 (i.e. 120.120 kHz)
    disp("Checking constants: the 1H-1H 1 angstrom dipolar coupling is "+prot_1a_hz+" Hz")
    
    f_1a = (mu_0 * hbar * gamma^2)/(4*pi*ang^3); % Calculate 1 angstrom F-F dipolar coupling (in radians!)
    
    dipcoup = f_1a;
            
    W = distance_matrix.^(-3)*dipcoup; % Homonuclear dipolar coupling strength

    % Neglect any coupling longer than the distance specified
    min_coupling = max_distance_computed^(-3)*dipcoup;
    W(W<min_coupling) = 0; 

    disp(min_coupling)
    
    
    Wsqu = W.^2; % Coupling squared
    
    for i = 1:N % Detailed balance
        Wsqu(i,i) = 0; % Zero the diagonal of the coupling matrix
        Wsqu_sums = sum(Wsqu,1); % Sum each column
        Wsqu(i,i) = -Wsqu_sums(i); % Set each diagonal to the negative sum of the column
    end   
            
    % Calculate the exchange matrix K
    K = 0.5*pi*Wsqu*powd*F0/1000000;  % Ignore (1-3*costheta^2) by replacing it with 0.8.
    propegator = expm((time_step/1000)*K); % Calculate the propegator
    
    current_matrix = eye(N); % Initial state is identity matrix

    all_matrices = zeros(N,N,num_time_points);

    for time_index = 1:num_time_points
        all_matrices(:,:,time_index) = current_matrix;
        current_matrix = propegator * current_matrix;
    end

    % Now, the magnetization that's detected in CODEX is whatever is left
    % on the diagonal
    out_magnetization = zeros(num_time_points,1);
    for time_index = 1:num_time_points
        out_magnetization(time_index) = trace(all_matrices(:,:,time_index))/N;
    end

    curve_out = out_magnetization;
end