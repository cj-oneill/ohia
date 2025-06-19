clear;clc; close all
%{
1. Define Radar Equation
    a) Model multiple different radar sats
        i) Sentinel 1
        ii) Nisar - 12 day repeat, S-band (10-12cm), 7m Azimuth range and
        3-24m slant range resolution, with indience angle of 33-47 deg. Open data
        access. 
        https://www.eoportal.org/satellite-missions/nisar#sensor-complement
        https://nisar.jpl.nasa.gov/mission/quick-facts/
        iii) ALOS-4 (PALSAR-3) - L - band (23cm), 3m Rg and 1m Az resolution, 8-70deg
        incidence angle range
        https://www.eorc.jaxa.jp/ALOS/en/alos-4/a4_sensor_e.htm
2. Send Pulse
    a) At fixed height, fixed speed, vary trees below
    b) Have to model velocity because you will increase aperature (SAR)
        i) Use SAR equations 
        ii) Model at angle of satellite
    c) Monte Carlo
        i) Take density of Ohia forest + generate random points
        ii) Variable percentage of healthy vs. dead trees
        iii) Model return using paper
3. Model Spatial + Temporal Resolution
    a) Take in S/C orbital parameters
    b) Frequency going over hawaii
    c) Period in which all of Hawaii is covered

Results: 
    1. Plots of example trees and covereage.
    2. Monte Carlo Results
    3. Monte Carlo of Percent Difference
    4. Table to evaulate simply

Assumptions: 
    1. Trees are circles
    2. Uniform distribution across x and y
    3. Trees radius cannot overlap on top of tree center
    4. RCS values stack (CHANGE this)
    5. Each resolution cell is average of RCS inside
    6. Best possible resolution (spotlight usually)

%colorbar = colormap(hot(N))
set(gca,'ColorOrder',cmap)

%}
%% Parameters

% Model Parameters
sample_distance = .1; % [m] distance between each sample. 0 would be most accurate.
rng(1); % sets seed for all rand() calls

%% Input Parameters
% Ohia Data
density_ohia = .01; % # [trees/m^2] GET GOOD NUMBER FOR HERE
percent_infected = .1; % [%] of OHIA infected GET GOOD NUMBER
ohia_diameters = [4.5 9]; % [m] average diameters of ohia crown
RCS_infected = 2;
RCS_healthy = 8;
RCS_base = 8; %RCS if radar hits ground

% Earth Data
R_E = 6378; %km
G = 6.67*10^-11;
M_E = 6*10^24; %kg

% % NISAR Given Data
% Altitude = 747; %km
% slant_res = 3; %m
% res_at = 7; % m along track
% angle_inc_max = 47; %degree
% swath = 150; %km

% NISAR Calculated Parameters
% res_ct = slant_res/sind(angle_inc_max); %cross track resolution
% V_circ = sqrt(G*M_E/(R_E*1000+Altitude*1000))/1000; %Velocity of NISAR


% Capella X-Band SAR Data
% res_at = .5; % m spotlight mode 5x5km
% slant_res = .3;% m spotlight mode
res_at = 1.2; % m stripmap mode 100x10km
slant_res = .75; % m stripmap mode
angle_inc_max = 50; % deg

% % Capella X-Band SAR Calculated Data
res_ct = slant_res/sind(angle_inc_max); %cross track resolution

%% Generating Ohia Sample
%{
1. Set limit of "Window"
2. Find Ohia density
3. Model returns
%}

n_ct_cells = 100; % number of resolution cells in cross-track
n_at_cells = 100; % number od resolution cells in along-track
n_cells = n_ct_cells*n_at_cells;

% Defining Window for Analysis
x_window_min = -res_ct*n_ct_cells/2;
x_window_max = res_ct*n_ct_cells/2;
x_window = [x_window_min x_window_max]; %cross-track window size [m]

y_window_min = -res_at*n_at_cells/2;
y_window_max = res_at*n_at_cells/2;
y_window = [y_window_min y_window_max]; %along track window size [m]

window_area = res_ct*n_ct_cells * res_at*n_ct_cells;

% Row for each tree, x and y location in column
N_trees = round(density_ohia*window_area); % [#] = [#/m^2 * m^2]
pt_tree = zeros(N_trees,3); %col 1/2: x/y loc; col 3: 0=healthy & 1=infected
pt_tree(:,1) = x_window_min + (x_window_max-x_window_min).*rand(N_trees,1); %x location of trees [m]
pt_tree(:,2) = y_window_min + (y_window_max-y_window_min).*rand(N_trees,1); %y location of trees [m]
pt_tree(:,3) = rand(N_trees,1) < percent_infected; % 1 = infected, 0 = healthy
pt_tree(:,4) = (ohia_diameters(1) + (ohia_diameters(2)-ohia_diameters(1)).*rand(N_trees,1))/2; %radius of tree

% Make sure that crown of one tree does not extend over trunk of another
for i = 1:N_trees
    valid = false;
    while ~valid
        % Generate a new location
        new_x = x_window_min + (x_window_max - x_window_min) * rand;
        new_y = y_window_min + (y_window_max - y_window_min) * rand;
        new_r = (ohia_diameters(1) + (ohia_diameters(2) - ohia_diameters(1)) * rand) / 2;
        % Check distance from all previous trees
        if i == 1
            % Always valid for the first tree
            valid = true;
        else
            existing_pts = pt_tree(1:i-1,1:2);
            dists = sqrt(sum((existing_pts - [new_x, new_y]).^2, 2));
            % Make sure all existing centers are outside the new radius
            if all(dists > new_r)
                valid = true;
            end
        end
    end
    % Assign the valid values
    pt_tree(i,1) = new_x;
    pt_tree(i,2) = new_y;
    pt_tree(i,4) = new_r;
end

%% RCS Calculation

% Start with creating resolution cells
cell_boundaries_ct = zeros(n_ct_cells+1, 1); % +1 because 2 cells have 3 lines
cell_boundaries_at = zeros(n_at_cells+1, 1); %

cell_ct_numbers = 0:n_ct_cells; % cross track cell
cell_at_numbers = 0:n_at_cells; % along track cell
[X_cell, Y_cell] = meshgrid(cell_ct_numbers, cell_ct_numbers);
xy_cell_pairs = [X_cell(:), Y_cell(:)];
% Cell grid each with associated average RCS
% [cell_x_number, cell_y_number, avg RCS in cell]

% cell_RCS = [X_cell(:), Y_cell(:), RCS_base*ones(size(xy_cell_pairs,1), 1)];
cell_RCS = RCS_base*ones(n_ct_cells,n_at_cells);

for i = 1:n_ct_cells+1
    % 1st boundary start at x_min goes to x_max
    cell_boundaries_ct(i) = x_window_min + res_ct * (i-1);
end
for i = 1:n_at_cells+1
    % 1st boundary start at x_min goes to x_max
    cell_boundaries_at(i) = y_window_min + res_at * (i-1);
end

% Take samples of RCS in each cell and average to find cell RCS value
% Matrix will be: (x_location, y_location, RCS of sample)
samples_ct_cell = 0:sample_distance:res_ct; % cross track cell samples
samples_at_cell = 0:sample_distance:res_at; % along track cell samples
n_samples_cell = length(samples_ct_cell)*length(samples_at_cell);

% Get all (x, y) combinations using meshgrid
[X, Y] = meshgrid(samples_ct_cell, samples_at_cell);
xy_pairs = [X(:)+x_window_min, Y(:)+y_window_min]; %start at bottom left
% Sample point for every place in a cell, with an associated RCS value for each 
% [x_location, y_location, RCS]
samples_cell = [xy_pairs, RCS_base*ones(size(xy_pairs,1), 1)];

% Going through each tree for each sample for each bin
for i = 1:n_ct_cells 
    for j = 1:n_at_cells
        % This is in each cell
        samples_cell(:,3) = RCS_base;
        for l = 1:n_samples_cell
            % This is each sample
            % i_sample = l + n_samples_cell*(l);
            x_sample = samples_cell(l,1);
            y_sample = samples_cell(l,2);
            RCS = samples_cell(l,3);
            for m = 1:N_trees
                % This is a sample looking at a tree
                x_tree = pt_tree(m,1);
                y_tree = pt_tree(m,2);
                is_infected = pt_tree(m,3);
                r_tree = pt_tree(m,4);

                %calculate distance between sample point and tree
                d = sqrt( (x_sample-x_tree)^2 + (y_sample-y_tree)^2 );

                %check if sample is under tree crown
                if d<r_tree
                    %Check if tree is healthy or infected
                    if is_infected == 1
                        RCS = RCS + RCS_infected;
                    else % healthy tree
                        RCS = RCS + RCS_healthy;
                    end
                end
            end
            samples_cell(l,3) = RCS; % To store in variable permanently
            % Leaving each sample
        end


        cell_RCS(i,j) = mean(samples_cell(:,3));
        % Leaving the bin
        samples_cell(:,2) = samples_cell(:,2) + res_at; 

        % formatSpec = 'j = %1.1f Moving AT: 1st sample x: %1.1f, y: %1.1f\n';
        % fprintf(formatSpec,j,samples_cell(1,1),samples_cell(1,2))
        % formatSpec = 'Moving AT: Last sample x: %1.1f, y: %1.1f\n\n';
        % fprintf(formatSpec,samples_cell(end,1),samples_cell(end,2))
    end
    %WORKING ON THIS
    samples_cell(:,2) = samples_cell(:,2) - n_at_cells*res_at; 
    % move samples over to next cross track cell
    samples_cell(:,1) = samples_cell(:,1) + res_ct;

    % formatSpec = 'Moving CT: Last sample x: %1.1f, y: %1.1f\n';
    % fprintf(formatSpec,samples_cell(1,1),samples_cell(1,2))
    % formatSpec = 'Moving CT: Last sample x: %1.1f, y: %1.1f\n\n';
    % fprintf(formatSpec,samples_cell(end,1),samples_cell(end,2))
end
cell_RCS_matched = flipud(cell_RCS'); %matched to plot

%% Sim Plot
hold on 
for i = 1:N_trees
    theta = linspace(0, 2*pi, 100);
    radius = pt_tree(i,4);
    xc = radius * cos(theta) + pt_tree(i,1); % x_location center of tree
    yc = radius * sin(theta) + pt_tree(i,2); % y_location center of tree
    if pt_tree(i,3) == 1 % if infected
        fill(xc, yc,[0.804, 0.522, 0.247]); % Brown Fill
    else % healthy
        fill(xc, yc,[0.133, 0.545, 0.133]); % Green Fill
    end
end
for i = 1:length(cell_boundaries_ct)
    xline(cell_boundaries_ct(i),'b',LineWidth=1)
end
for i = 1:length(cell_boundaries_at)
    yline(cell_boundaries_at(i),'b',LineWidth=1)
end

xlim(x_window);
ylim(y_window)
title('Randomnly Distributed Ohia')
xlabel('Cross-Track [m]')
ylabel('Along-Track [m]')
hold off

%% Image Plot
%% Plot Rectangles Colored by cell_RCS (Custom Green-Brown Color Map)

figure;
hold on

% Custom colormap from dark brown to bright green
Ncolors = 256;
brown = [0.396, 0.263, 0.129]; % dark brown
green = [0.133, 0.545, 0.133]; % bright green
cmap_custom = [linspace(brown(1), green(1), Ncolors)', ...
               linspace(brown(2), green(2), Ncolors)', ...
               linspace(brown(3), green(3), Ncolors)'];

% Normalize RCS values
min_RCS = min(cell_RCS(:));
max_RCS = max(cell_RCS(:));
norm_RCS = (cell_RCS - min_RCS) / (max_RCS - min_RCS); % Normalize to [0,1]

% Draw each resolution cell as a colored rectangle
for i = 1:n_ct_cells
    for j = 1:n_at_cells
        % Cell corners
        x_left = x_window_min + (i-1)*res_ct;
        x_right = x_left + res_ct;
        y_bottom = y_window_min + (j-1)*res_at;
        y_top = y_bottom + res_at;

        % Get corresponding color from custom colormap
        color_idx = round(norm_RCS(i,j) * (Ncolors-1)) + 1;
        rect_color = cmap_custom(color_idx, :);

        % Draw rectangle
        x_rect = [x_left x_right x_right x_left];
        y_rect = [y_bottom y_bottom y_top y_top];
        fill(x_rect, y_rect, rect_color, 'EdgeColor', 'none');
    end
end

xlim(x_window);
ylim(y_window);
axis equal
colormap(cmap_custom)
cb = colorbar;
cb.Label.String = 'Average RCS';
clim([min_RCS, max_RCS])

title('SAR Resolution Cells Colored by Average RCS')
xlabel('Cross-Track [m]')
ylabel('Along-Track [m]')
hold off

