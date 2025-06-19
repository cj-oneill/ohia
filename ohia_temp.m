clear; clc; close all

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
    7. No diurnal water content changes... See the dutch masters thesis...
    8. Negligible atmospheric impacts


%colorbar = colormap(hot(N))
set(gca,'ColorOrder',cmap)
%}

%% Parameters
sample_distance = .1; % Distance between samples inside each resolution cell [m]
rng(); % Seed the random number generator for reproducibility

%% Input Parameters
density_ohia = .08; % Tree density [trees/m^2]
percent_infected = .1; % Percentage of infected trees
ohia_diameters = [4.5 9]; % Min and max tree diameter [m]
RCS_infected = 2; % RCS contribution from infected trees
RCS_healthy = 8; % RCS contribution from healthy trees
RCS_base = 8; % Base RCS in each cell

% Radar Parameters
res_at = 1.2; % Along-track resolution [m]
slant_res = .75; % Slant-range resolution [m]
angle_inc_max = 50; % Maximum incidence angle [deg]
res_ct = slant_res/sind(angle_inc_max); % Cross-track resolution

%% Generate Window & Ohia Tree Distribution
n_ct_cells = 100; % Number of cross-track cells
n_at_cells = 100; % Number of along-track cells

% Define image window dimensions
x_window_min = -res_ct*n_ct_cells/2;
x_window_max = res_ct*n_ct_cells/2;
y_window_min = -res_at*n_at_cells/2;
y_window_max = res_at*n_at_cells/2;
x_window = [x_window_min x_window_max];
y_window = [y_window_min y_window_max];
window_area = res_ct*n_ct_cells * res_at*n_at_cells;

% Place trees randomly
N_trees = round(density_ohia*window_area);
pt_tree = zeros(N_trees,4); % Columns: x, y, infected flag, radius
pt_tree(:,1) = x_window_min + (x_window_max-x_window_min).*rand(N_trees,1);
pt_tree(:,2) = y_window_min + (y_window_max-y_window_min).*rand(N_trees,1);
pt_tree(:,3) = rand(N_trees,1) < percent_infected;
pt_tree(:,4) = (ohia_diameters(1) + (ohia_diameters(2)-ohia_diameters(1)).*rand(N_trees,1))/2;

% Reposition trees to avoid overlap
for i = 1:N_trees
    valid = false;
    while ~valid
        new_x = x_window_min + (x_window_max - x_window_min) * rand;
        new_y = y_window_min + (y_window_max - y_window_min) * rand;
        new_r = (ohia_diameters(1) + (ohia_diameters(2) - ohia_diameters(1)) * rand) / 2;
        if i == 1
            valid = true;
        else
            existing_pts = pt_tree(1:i-1,1:2);
            dists = sqrt(sum((existing_pts - [new_x, new_y]).^2, 2));
            if all(dists > new_r)
                valid = true;
            end
        end
    end
    pt_tree(i,1:2) = [new_x, new_y];
    pt_tree(i,4) = new_r;
end

%% RCS Calculation (Cell Sampling)
cell_boundaries_ct = linspace(x_window_min, x_window_max, n_ct_cells+1);
cell_boundaries_at = linspace(y_window_min, y_window_max, n_at_cells+1);
cell_RCS = RCS_base*ones(n_ct_cells,n_at_cells); % Preallocate RCS grid

% Create samples inside each cell
samples_ct_cell = 0:sample_distance:res_ct;
samples_at_cell = 0:sample_distance:res_at;
[X, Y] = meshgrid(samples_ct_cell, samples_at_cell);
xy_pairs = [X(:)+x_window_min, Y(:)+y_window_min];
samples_cell = [xy_pairs, RCS_base*ones(size(xy_pairs,1), 1)];
n_samples_cell = size(samples_cell, 1);

% Loop through resolution cells
for i = 1:n_ct_cells 
    for j = 1:n_at_cells
        samples_cell(:,3) = RCS_base; % Reset RCS values
        for l = 1:n_samples_cell
            x_sample = samples_cell(l,1);
            y_sample = samples_cell(l,2);
            RCS = samples_cell(l,3);
            for m = 1:N_trees
                x_tree = pt_tree(m,1);
                y_tree = pt_tree(m,2);
                is_infected = pt_tree(m,3);
                r_tree = pt_tree(m,4);
                d = sqrt((x_sample-x_tree)^2 + (y_sample-y_tree)^2);
                if d < r_tree
                    if is_infected == 1
                        RCS = RCS + RCS_infected;
                    else
                        RCS = RCS + RCS_healthy;
                    end
                end
            end
            samples_cell(l,3) = RCS;
        end
        cell_RCS(i,j) = mean(samples_cell(:,3)); % Average RCS in cell
        samples_cell(:,2) = samples_cell(:,2) + res_at; % Shift Y
    end
    samples_cell(:,2) = samples_cell(:,2) - n_at_cells*res_at;
    samples_cell(:,1) = samples_cell(:,1) + res_ct; % Shift X
end

%% Tree + Grid Plot
figure;
hold on
for i = 1:N_trees
    theta = linspace(0, 2*pi, 100);
    radius = pt_tree(i,4);
    xc = radius * cos(theta) + pt_tree(i,1);
    yc = radius * sin(theta) + pt_tree(i,2);
    if pt_tree(i,3) == 1
        fill(xc, yc, [0.804, 0.522, 0.247]); % brown for infected
    else
        fill(xc, yc, [0.133, 0.545, 0.133]); % green for healthy
    end
end
% Draw resolution grid
for i = 1:length(cell_boundaries_ct)
    xline(cell_boundaries_ct(i), 'b', 'LineWidth', 1)
end
for i = 1:length(cell_boundaries_at)
    yline(cell_boundaries_at(i), 'b', 'LineWidth', 1)
end
axis tight
axis equal
set(gca, 'XLim', x_window, 'YLim', y_window)
set(gcf, 'Position', [80, 80, 600, 600])
xlabel('Cross-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Along-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Simulated Ohia Forest with Resolution Cells', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'FontName', 'Arial');
hold off

%% RCS Heatmap Plot (Rectangles + Colorbar)
figure;
hold on

% Custom colormap: dark brown to bright green
Ncolors = 256;
% color_1 = [0.396, 0.263, 0.129]; % dark brown
% color_2 = [0.133, 0.545, 0.133]; % bright green

color_1 = [0, 0, 0]; % black
%% 
color_2 = [1, 1, 1]; % white

cmap_custom = [linspace(color_1(1), color_2(1), Ncolors)', ...
               linspace(color_1(2), color_2(2), Ncolors)', ...
               linspace(color_1(3), color_2(3), Ncolors)'];

% Normalize RCS values for color mapping
min_RCS = min(cell_RCS(:));
max_RCS = max(cell_RCS(:));
norm_RCS = (cell_RCS - min_RCS) / (max_RCS - min_RCS);

% Draw colored resolution cells
for i = 1:n_ct_cells
    for j = 1:n_at_cells
        x_left = x_window_min + (i-1)*res_ct;
        x_right = x_left + res_ct;
        y_bottom = y_window_min + (j-1)*res_at;
        y_top = y_bottom + res_at;
        color_idx = round(norm_RCS(i,j) * (Ncolors-1)) + 1;
        rect_color = cmap_custom(color_idx, :);
        fill([x_left x_right x_right x_left], [y_bottom y_bottom y_top y_top], rect_color, 'EdgeColor', 'none');
    end
end

colormap(cmap_custom)
cb = colorbar;
cb.Label.String = 'Average RCS';  
cb.Label.FontSize = 12;  
cb.Label.FontWeight = 'bold';  
cb.Label.Interpreter = 'latex';
clim([min_RCS, max_RCS])
axis equal
axis([x_window y_window])
set(gcf, 'Position', [80, 80, 600, 600])
xlabel('Cross-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Along-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Simulated Ohia Radar Coverage', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
set(gca, 'FontSize', 12, 'FontName', 'Arial');
hold off
