%-------------------------------------------------------------------------%
% ClusteringMonteCarlo.m
% Simulated clustering of oligomers given a pairwise interaction potential
% Westley W. Wu, Noah H. Somberg, Mei Hong
% Written in MATLAB R2021b
% July 2022
% This script takes a specified number of pores and places them randomly on
% a 1000 A x 1000 A area. It then applies the interaction potential and
% simulates with a Metropolis Monte Carlo method. The interaction potential
% is read from the file 'interaction_potential.txt'. This simulation uses
% periodic boundary conditions and outputs the final coordinates of the
% center of each oligomer in 'finalCenters.csv'
%-------------------------------------------------------------------------%

%%%%%%%% Begin main program %%%%%%%%

%%% Initialize a random list of ordered pairs to serve as an initial random
%%% distribution of pores

% How many pores we place — for 1:17 P:L ratio in a 1000 A x 1000 A box,
% 349 pores are appropriate
npores = 326; 
% Generates a membrane patch with 1:17 P:L ratio given the number of pores
sl = sqrt(1000000*(npores/349));
% Initial pore center locations are created. MATLAB randomly picks these
% center coordinates based on the number of pores specified by npores. Note
% that oligomer overlap is allowed
coords = sl*rand(npores,2);
% Saves a record of these initial centers just for future reference
writematrix(coords, 'originalCenters.csv');

% Number of iterations of the Metropolis Monte Carlo algorithm we want to
% simulate. In other words, how many individual oligomer positional changes
% we consider making
ntimes = 15000;

samp = 50; % Take snapshots every 50 iters
snapshots = zeros(326,2,floor(ntimes/samp));
snapshots(:,:,1) = coords;
iSnap = 2;
%%% For later plotting of average nearest neighbor (NN) distances

% Average nearest neighbor (NN) distances table: first column will be the
% iteration the simulation is on. The second column will be the NN distance
% averaged over all npores oligomers.
nnVtime = zeros(ntimes+1, 2);
% Computes an average NN distance for the initial pore setup prior to any
% Monte Carlo simulations
nnVtime(1, :) = [0 ANNDistance(coords, sl)];

% % For each cycle, we attempt to move a randomly-chosen oligomer to a randomly-
% % chosen position. We then use the change in energy between the current and
% % proposed configurations to decide if we should accept or reject the proposed
% % configuration. Periodic boundary conditions are used.
for i = 1:ntimes
    currentRow = randi(npores); % Chooses a random oligomer to try and move

    % coordsFinal always refers to the proposed state of the system. Here,
    % we reset the proposed state of the system to the most recently-accepted
    % state at the start of each loop. When proposing a new state later on
    % in this code, we use the current state (=coordsFinal at this point) as
    % a foundation onto which we make a small modification to arrive at our
    % new proposed state
    coordsFinal = coords;
    
    %%% Gets positional information about randomly-chosen point
    current_x = coords(currentRow,1);
    current_y = coords(currentRow,2);

    %%% Compute the energy of our initial state

    % Finds nearest centers relevant enough to include in an energy calculation
    % of the current state. Centers that are not considered relevant are
    % those whose distances are too far away from the pore we plan to move
    % (and hence contribute 0 pairwise interaction energy)
    nearOriginal = nearMe(current_x, current_y, coords, sl);
    % Actually computes the energy between chosen center and its close neighbors
    E_i = computeEnergy(coords(currentRow,:), nearOriginal); 
    
    %%% Picks up the pore we randomly chose before and moves it to a random
    %%% destination to form the proposed state

    % Randomly chooses a proposed destination
    finalCoordsForPoint = sl*rand(1,2);
    % Changes the position of the chosen pore only in coordsFinal (note that
    % coords, the current state, remains unchanged)
    coordsFinal(currentRow,:) = finalCoordsForPoint;
    % Positional information about the point's new location
    current_x_f = coordsFinal(currentRow,1);
    current_y_f = coordsFinal(currentRow,2);

    %%% Compute the energy of our proposed state, as we've done before with
    %%% the current state

    % Finds nearest centers relevant enough to include in an energy calculation
    % of the proposed state. Centers that are not considered relevant are
    % those whose distances are too far away from the pore we plan to move
    % (and hence contribute 0 pairwise interaction energy)
    nearFinal = nearMe(current_x_f, current_y_f, coordsFinal, sl);
    % Actually computes the energy between chosen center and its close neighbors
    E_f = computeEnergy(coordsFinal(currentRow,:), nearFinal);
   
    % Given the energies of the current and propsed states, we make a
    % decision: do we accept the new state? By the detailed balance conditon,
    % if the energy of the proposed state is higher (less favorable) than
    % the current state, the proposed state is accepted with probability
    % e^(E_i-E_f). If the energy of the proposed state is lower (more
    % favorable) than the current state, the proposed state is always
    % accepted.
    if E_i-E_f < 0 % Energy of proposed is less favorable
        if rand <= exp(E_i-E_f) % Acceptance probability implementation
            % If accepted, change the current state into the proposed state
            coords = coordsFinal;
        end
        % If not accepted, do nothing on this iteration
    else % Energy of proposed is more favorable
        % Change the current state to reflect the new proposed state
        coords = coordsFinal;
    end

    if mod(i,samp) == 0
        snapshots(:,:,iSnap) = coords;
        iSnap = iSnap + 1;
    end

    % Updates the average NN computation datatable with first entry =
    % iteration number and second entry = average NN distance
    nnVtime(i+1, :) = [i ANNDistance(coords, sl)];

    % Prints which iteration we are on, just for reference and can delete
    % this line if so desired
    if mod(i,100) == 0
        disp(i);
    end
end

% To display periodic boundary conditions, we take our current box and
% create a 3x3 grid with it (only for figure 4 below)
q1 = [coords(:,1)-sl coords(:,2)+sl];
q2 = [coords(:,1) coords(:,2)+sl];
q3 = [coords(:,1)+sl coords(:,2)+sl];
q4 = [coords(:,1)-sl coords(:,2)];
q5 = coords;
q6 = [coords(:,1)+sl coords(:,2)];
q7 = [coords(:,1)-sl coords(:,2)-sl];
q8 = [coords(:,1) coords(:,2)-sl];
q9 = [coords(:,1)+sl coords(:,2)-sl];
expandedGrid = [q1;q2;q3;q4;q5;q6;q7;q8;q9];

%%% Saves results

% Final oligomer center positions with periodic boundary conditions
writematrix(expandedGrid, 'expandedFinalCenters.csv');
% Final oligomer centers without a display of periodicity
writematrix(coords, 'finalCenters.csv');
% Monte Carlo algorithm iteration vs average NN datatable
writematrix(nnVtime, 'avg_NN_dist_v_time.csv'); 

%%% Figures below, for descriptions of each see their titles

figure(4);
plot(nnVtime(:,1),nnVtime(:,2));
title('Average nearest neighbor distance over time');
xlabel('Timestep');
ylabel('Mean NN distance (angstrom)');
xlim([0 ntimes]);

figure(5);
polygonPoints(5, 8.8, 0);
xlim([0 sl]);
ylim([0 sl]);
title('Final centers');

figure(6);
polygonPoints(5, 8, 1);
xlim([0 sl]);
ylim([0 sl]);
title('Original centers');

figure(7);
polygonPoints(5, 8, 2);
xlim([-sl 2*sl]);
ylim([-sl 2*sl]);
title('Final centers w/ periodic boundary conditions');

figure(8);
x = movmean(nnVtime(:,1),500);
y = movmean(nnVtime(:,2),500);
plot(x,y)
title('Average nearest neighbor distance over time');
xlabel('Timestep');
ylabel('Moving average of mean NN distance (angstrom)');


%%%%%%%% Functions used in main simulation begin here %%%%%%%%

% Computes the total energy of all pairwise interactions between a given
% oligomer center's coordinates ('pt') and every other oligomer center
% coordinate which appears in a list of pt's nearest neighbors ('nearby')
function energy = computeEnergy(pt, nearby)
    % Vectorized computation of distance from pt to every other point in
    % nearby; dist ends up being a list of all pairwise distances
    pointMatrix = repmat(pt, size(nearby,1), 1);
    dist = sqrt(sum(((pointMatrix-nearby).^2),2));

    % Computes the total energy from the list of pairwise distances
    if size(dist,1)>=1 % Verifies that pt actually has neighbors
        dist = double(dist); % Converts all distances to double just to be safe
        % Runs distance list through potential energy function, end up with
        % a list of all pairwise energies
        energies = potential(dist);
        % Adds together all the individual pairwise energies computed in the line above
        energy = sum(energies,'all');
    else
        % If pt doesn't have neighbors, there is no total pairwise energy.
        % So, we save effort and don't compute it
        energy = 0.0;
    end
end

% Potential function, r is center-center distance in angstroms; function
% output is the energy in units of kB*T at a center-center distance r
function energy = potential(r)
    values = readtable('interaction_potential.txt'); % Reads potential from file
    % First column of text file contains distances in angstroms
    d = table2array(values(:,1));
    % Second column of text file contains energies in units of kB*T
    y = table2array(values(:,2));
    % Datafile is a discrete table, so we use MATLAB's modified Akima method
    % to interpolate between points; pp ends up being the pairwise potential
    % energy function we want to use
    pp = griddedInterpolant(d,y,'makima');
    % Plugs in distance r into function pp to compute a numerical value of
    % the energy
    energy = pp(double(r));
end

% Finds the center coordinates of all oligomers at distances close enough to
% a particular pore center (x,y) such that the pairwise interaction potential
% with (x,y) is possibly nonzero. The side length of non-periodic membrane
% patch ('sl') is needed to implement periodic boundary conditions.
% This function ensures we don't compute more energies than we have to,
% since beyond center-center 90 A, all pairwise energies become 0 anyway.
% This function exists simply to improve simulation speed
function expandedGrid = nearMe(x,y,coords,sl)
    
    % Initiates 9 copies of the grid in 3x3 arrangement to implement periodic
    % boundary conditions; expandedGrid is a list of all pore centers in
    % this 3x3 repeated grid
    q1 = [coords(:,1)-sl coords(:,2)+sl];
    q2 = [coords(:,1) coords(:,2)+sl];
    q3 = [coords(:,1)+sl coords(:,2)+sl];
    q4 = [coords(:,1)-sl coords(:,2)];
    q5 = coords;
    q6 = [coords(:,1)+sl coords(:,2)];
    q7 = [coords(:,1)-sl coords(:,2)-sl];
    q8 = [coords(:,1) coords(:,2)-sl];
    q9 = [coords(:,1)+sl coords(:,2)-sl];
    expandedGrid = [q1;q2;q3;q4;q5;q6;q7;q8;q9];
    
    % Keeps track of which centers from expandedGrid we want to remove from
    % consideration. By default, remove the center unless told otherwise
    rowsToRemove = true(size(expandedGrid,1), 1);

    % Keeps track of how many points in expandedGrid have exact coorindates
    % as our query point (x,y); this is only applicable to the extremely
    % unlikely edge case of different pores having centers with identical
    % coordinates. Even though this is a nearly-impossible edge case, we
    % still want the program to handle it just in case
    same_counter = 0; 

    % Removes all centers in expandedGrid satisfying either of the following
    % conditions:
    % (1) the center is exactly (x,y), or
    % (2) located beyond a 200 A x 200 A square centered on (x,y).
    % Note that condition (1) is present to avoid accidentally computing the
    % nonexistent pairwise potential between a pore with itself and that (2)
    % is more than sufficient to capture all points within a 90 A radius of
    % (x,y)
    for row = 1:size(expandedGrid,1) % Cycle through all expandedGrid centers
        if (expandedGrid(row,1) == x) && (expandedGrid(row,2) == y) % If the coordinate is the same as the query point...
            same_counter = same_counter+1; % Keep track of how many there are
            % If the same_counter = 1, this is to be expected, meaning that
            % in our expandedGrid, there is exactly one point with center
            % coordinates (x,y), namely (x,y) itself. If same_counter > 1,
            % we know that there are same_counter – 1 other pores sharing the
            % same center (x,y) as our query pore
        end
        % Tests conditions (1) and (2) at the same time
        if (expandedGrid(row,1) >= x-100) && (expandedGrid(row,1) <= x+100)...
                && (expandedGrid(row,2) >= y-100) && (expandedGrid(row,2) <= y+100)...
                && (expandedGrid(row,1) ~= x) && (expandedGrid(row,2) ~= y)
            % Do not remove center if it lies within the 200 A x 200 A
            % square of interest AND the center is not (x,y) itself
            rowsToRemove(row,1) = false;
        end
    end

    expandedGrid(rowsToRemove, :) = []; % Deletes all centers to be removed

    % Accounts for the unlikely edge case of multiple different pores
    % sharing the same center as (x, y)
    temp = repmat([x y], same_counter-1, 1);
    expandedGrid = [expandedGrid; temp];

    % We end up with a list of all pore centers within a 200 A x 200 A box
    % of the query point in expandedGrid, EXCEPT for the coordinates
    % corresponding to the query point itself
end

% Plots randomly-angled symmetric n-gons from centers files generated
% from the Monte Carlo simulation. Which center file gets plotted is dictated
% by the input parameter 'p.' The code here is ONLY for purposes of
% visualization after the simulation is complete. The polygon orientations,
% vertices, and distance matrix used in CODEX are generated in our CODEX
% pentagon-plotting script, not here. The output to the entire Monte Carlo
% simulation is only a file of oligomer center coordinates, not a
% distance matrix of n-gon vertices
function polygonPoints(n,s,p) % n = oligomer number; s = side length; p = which file to plot
    
    if  p == 0 % Plots final centers
        centers = readmatrix("finalCenters.csv"); 
    elseif p == 1 % Plots original centers before any iterations of Monte Carlo
        centers = readmatrix("originalCenters.csv"); 
    else % Plots final centers zoomed out showing periodicity
        centers = readmatrix("expandedFinalCenters.csv");
    end

    % randVerts to contain the coordinates of the vertices of a regular n-gon
    % with center at the origin, (0,0)
    randVerts = zeros(1,n);
    circumradius = (double(s)/2)/(sin(pi/n)); % Computes circumradius

    % Cycles through pore centers, generating a randomly-angled regular n-gon
    % centered on each specified center-coordinate
    for j = 1:size(centers,1) 
        % Make a randomly-oriented n-gon centered at the origin...
        randStart = rand*2*pi/n;
        for i = 1:n
            randVerts(i) = randStart+2*pi*i/n;
        end

        % ... and shift the n-gon made above to where it is supposed to be
        % centered
        x_coord = centers(j,1)+circumradius*cos(randVerts);
        y_coord = centers(j,2)+circumradius*sin(randVerts);

        % Plots the entire grid of pentagons
        patch(x_coord, y_coord, 'blue');
    end
end

% Computes an average NN distance for a given list of center coordinates
% ('coords'). The side length of non-periodic membrane patch ('sl') is needed
% to implement periodic boundary conditions
function dist = ANNDistance(coords, sl)
   
    % Initiates 9 copies of the grid in 3x3 arrangement to implement periodic
    % boundary conditions; expandedGrid is a list of all pore centers in
    % this 3x3 repeated grid
    q1 = [coords(:,1)-sl coords(:,2)+sl];
    q2 = [coords(:,1) coords(:,2)+sl];
    q3 = [coords(:,1)+sl coords(:,2)+sl];
    q4 = [coords(:,1)-sl coords(:,2)];
    q5 = coords;
    q6 = [coords(:,1)+sl coords(:,2)];
    q7 = [coords(:,1)-sl coords(:,2)-sl];
    q8 = [coords(:,1) coords(:,2)-sl];
    q9 = [coords(:,1)+sl coords(:,2)-sl];
    expandedGrid = [q1;q2;q3;q4;q5;q6;q7;q8;q9];

    % Running sum of NN distances
    dist = 0;
    % Finds the nearest neighbors to each coords point in expandedGrid
    idx = knnsearch(expandedGrid,coords,"K",2);
    for i = 1:size(coords,1) % Cycles through each pore in ONLY the non-expanded grid...
        % ... but NN computations take into account centers that appear as a
        % result of periodic boundary conditions. The expression under the
        % square root computes the NN distance for pore i
        dist = dist + sqrt((coords(i,1)-expandedGrid(idx(i,2),1))^2 +...
            (coords(i,2)-expandedGrid(idx(i,2),2))^2);
    end
    dist = dist/size(coords,1); % Averages running over total number of pores
end

%-------------------------------------------------------------------------%
% End of ClusteringMonteCarlo.m
%-------------------------------------------------------------------------%