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

%Sim Inputs
time = 1; % 0: no time sim, 1: time sim;
std = 0; % 0: no std sim, 1: std sim;
sample_distance = .1; % Distance between samples inside each resolution cell [m]
n_ct_cells = 30;
n_at_cells = 30;
rng(2); % Seed the random number generator for reproducibility

% Export Inputs
export_std = 0; % 0: no export, 1: export; for standard Sim/SAR export
export_time = 1; % for 3 col time dependence 
export_inf_chance = 1; % for infection chance versus distance
export_dielecric = 1; % for dielectric vs. time infected
export_refl_t = 0; % for reflectivity vs. time
folder = 'Figures';
fileName = 'Ohia';

% Ohia Data
density_ohia = .125; % # [trees/m^2] GET GOOD NUMBER FOR HERE
percent_infected = .1; % [%] of OHIA infected GET GOOD NUMBER
ohia_diameters = [4.5 9]; % [m] average diameters of ohia crown
dielectric_infected = 4;
dielectric_healthy = 30;
dielectric_std = 9; %dielectric 1 sigma error
time_to_die = 14; % days

RCS_base = 0; %RCS if radar hits ground

% Spread Modeling
max_infection_chance = .5; % in percentage 
dist_no_spread = 20; %m

% Time
t_0 = 0; %days
t_f = 24; %days (12 days per pass)
t_check_day = 1; % number of times simulation run per day?


% Earth Data
R_E = 6378; %km
G = 6.67*10^-11;
M_E = 6*10^24; %kg

% % Sentinel-1 Given Data
% res_ct = 5; %m
% res_at = 5; % m along track
% wavelength = .05; %S-Band


% % NISAR Calculated Parameters
% res_ct = slant_res/sind(angle_inc_max); %cross track resolution
% V_circ = sqrt(G*M_E/(R_E*1000+Altitude*1000))/1000; %Velocity of NISAR
% 
% % NISAR Given Data
% Altitude = 747; %km
% slant_res = 3; %m
% res_at = 7; % m along track
% angle_inc_max = 47; %degree
% swath = 150; %km
% wavelength = .24; %S-Band
% % wavelength = .1; %L-Band


% % NISAR Calculated Parameters
% res_ct = slant_res/sind(angle_inc_max); %cross track resolution
% V_circ = sqrt(G*M_E/(R_E*1000+Altitude*1000))/1000; %Velocity of NISAR


% Capella X-Band SAR Data
% res_at = .5; % m spotlight mode 5x5km
% slant_res = .3;% m spotlight mode
res_at = 1.2; % m stripmap mode 100x10km
slant_res = .75; % m stripmap mode
angle_inc_max = 50; % deg
wavelength = .03; %m CHECK

% % Capella X-Band SAR Calculated Data
res_ct = slant_res/sind(angle_inc_max); %cross track resolution

%% Time

if time == 1
    % Generate points with only 1% infected
    percent_infected = .03; 
    % Generating trees
    num_t_sim = (t_f-t_0) * t_check_day; %total # of days * checks per day
    %pt tree is the local simulation always changing
    %time_pt_tree is savings copies of pt_tree at various times
    [pt_tree,N_trees,x_window,y_window] = tree_gen(n_ct_cells,n_at_cells,res_ct,res_at,density_ohia,percent_infected,ohia_diameters);
    time_pt_tree = zeros(num_t_sim, N_trees,5); % Columns: x, y, infected flag, radius


    %Create cell Boundaries
    cell_boundaries_ct = linspace(x_window(1), x_window(2), n_ct_cells+1);
    cell_boundaries_at = linspace(y_window(1), y_window(2), n_at_cells+1);
    time_cell_RCS = RCS_base*ones(num_t_sim,n_ct_cells,n_at_cells);

    t = linspace(t_0,t_f,num_t_sim);

    for s = 1:num_t_sim % # of simulations
       if s == 1 % 1st simulation
           % dont do anything
           time_pt_tree(1,:,:) = pt_tree;
       else
            % Simulate spread of ROD
            for i = 1:N_trees
                % Use infection spread function to simualte spread based on
                % distance from any other tree
                if pt_tree(i,3) == 1 % if host tree is infected
                    x_host = pt_tree(i,1);
                    y_host = pt_tree(i,2);
                    if mod(s,t_check_day) == 0 % once per day
                        %increase # of days tree has been infected for
                        pt_tree(i,5) = pt_tree(i,5) + 1;
                    end
                    for j = 1:N_trees
                        if pt_tree (j, 3) == 0 % target tree is healthy
                            x_trgt = pt_tree(j,1);
                            y_trgt = pt_tree(j,1);
                            % infection chance is function of distance
                            d = sqrt((x_host-x_trgt)^2 + (y_host-y_trgt)^2);
                            chance = inf_sprd(d,pt_tree(i,3));

                            % assigns if trgt tree is infected
                            pt_tree(j,3) = rand(1) < chance;
                        end
                    end
                end
            end
       end
       time_pt_tree(s,:,:) = pt_tree;
       time_cell_RCS(s,:,:) = get_cell_RCS(x_window, y_window, n_ct_cells, n_at_cells, res_ct, res_at, sample_distance, RCS_base, pt_tree, dielectric_healthy, dielectric_infected, dielectric_std, wavelength,time_to_die);
    end


    %% Time Plot
    figure()
    t_to_plot = [1 12 24]; % days for columns of sims
    sim_i = t_to_plot*t_check_day;
    for i = 1:length(t_to_plot)
        % Plot the sim
        subplot(2,length(t_to_plot),i)
        hold on
        for j = 1:N_trees
            theta = linspace(0, 2*pi, 100);
            radius = time_pt_tree(sim_i(i),j,4);
            xc = radius * cos(theta) + time_pt_tree(sim_i(i),j,1);
            yc = radius * sin(theta) + time_pt_tree(sim_i(i),j,2);
            if time_pt_tree(sim_i(i),j,3) == 1 %If j tree in i sim is infected
                fill(xc, yc, [0.804, 0.522, 0.247]); % brown for infected
            else
                fill(xc, yc, [0.133, 0.545, 0.133]); % green for healthy
            end
        end
        % Draw resolution grid
        for j = 1:length(cell_boundaries_ct)
            xline(cell_boundaries_ct(j), 'b', 'LineWidth', 1)
        end
        for j = 1:length(cell_boundaries_at)
            yline(cell_boundaries_at(j), 'b', 'LineWidth', 1)
        end
        axis tight
        axis equal
        set(gca, 'XLim', x_window, 'YLim', y_window)
        % set(gcf, 'Position', [80, 80, 600, 600])
        xlabel('[m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
        ylabel('[m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        hold off
    

        % Plot the SAR image
        subplot(2,length(t_to_plot),i+length(t_to_plot)) %bottom row
        hold on

        % Custom colormap: dark brown to bright green
        Ncolors = 256;
        color_1 = [0, 0, 0]; % black
        color_2 = [1, 1, 1]; % white
        cmap_custom = [linspace(color_1(1), color_2(1), Ncolors)', ...
                       linspace(color_1(2), color_2(2), Ncolors)', ...
                       linspace(color_1(3), color_2(3), Ncolors)'];
        
        % Normalize RCS values for color mapping
        rcs_temp_sim = squeeze(time_cell_RCS(sim_i(i),:,:));
        min_RCS = min(rcs_temp_sim(:));
        max_RCS = max(rcs_temp_sim(:));
        
        %Check if want this
        norm_RCS = squeeze((time_cell_RCS(sim_i(i),:,:) - min_RCS) / (max_RCS - min_RCS));
        
        % Draw colored resolution cells
        for k = 1:n_ct_cells
            for j = 1:n_at_cells
                x_left = x_window(1) + (k-1)*res_ct;
                x_right = x_left + res_ct;
                y_bottom = y_window(1) + (j-1)*res_at;
                y_top = y_bottom + res_at;
                color_idx = round(norm_RCS(k,j) * (Ncolors-1)) + 1;
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
        xlabel('[m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
        ylabel('[m]', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
        set(gca, 'FontSize', 12, 'FontName', 'Arial');
        hold off
        
    end
    set(gcf,'Position',[500,500,1200,600]);

    if export_time == 1
        exportgraphics(gcf, [folder,'\', fileName,'time', '.pdf'],'ContentType','vector')
    end

end
% Generate Wirndow & Ohia Tree Distribution
% 
% [pt_tree,N_trees,x_window,y_window] = tree_gen(n_ct_cells,n_at_cells,res_ct,res_at,density_ohia,percent_infected,ohia_diameters);
% 
% 
% %% RCS Calculation (Cell Sampling)
% 
% %Calculate RCS of generatred trees
% cell_RCS = get_cell_RCS(x_window, y_window, n_ct_cells, n_at_cells, res_ct, res_at, sample_distance, RCS_base, pt_tree, dielectric_healthy, dielectric_infected, dielectric_std, wavelength);
% 
% 
% %% Tree + Grid Plot
% figure;
% hold on
% for i = 1:N_trees
%     theta = linspace(0, 2*pi, 100);
%     radius = pt_tree(i,4);
%     xc = radius * cos(theta) + pt_tree(i,1);
%     yc = radius * sin(theta) + pt_tree(i,2);
%     if pt_tree(i,3) == 1
%         fill(xc, yc, [0.804, 0.522, 0.247]); % brown for infected
%     else
%         fill(xc, yc, [0.133, 0.545, 0.133]); % green for healthy
%     end
% end
% % Draw resolution grid
% cell_boundaries_ct = linspace(x_window(1), x_window(2), n_ct_cells+1);
% cell_boundaries_at = linspace(y_window(1), y_window(2), n_at_cells+1);
% for i = 1:length(cell_boundaries_ct)
%     xline(cell_boundaries_ct(i), 'b', 'LineWidth', 1)
% end
% for i = 1:length(cell_boundaries_at)
%     yline(cell_boundaries_at(i), 'b', 'LineWidth', 1)
% end
% axis tight
% axis equal
% set(gca, 'XLim', x_window, 'YLim', y_window)
% set(gcf, 'Position', [80, 80, 600, 600])
% xlabel('Cross-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('Along-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% title('Simulated Ohia Forest with Resolution Cells', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
% set(gca, 'FontSize', 12, 'FontName', 'Arial');
% hold off
% 
% if export_std == 1
%     exportgraphics(gcf, [folder,'\',fileName,'.pdf'], 'ContentType', 'vector')
% end
% 
% 
% %% RCS Heatmap Plot (Rectangles + Colorbar)
% figure;
% hold on
% 
% % Custom colormap: dark brown to bright green
% Ncolors = 256;
% % color_1 = [0.396, 0.263, 0.129]; % dark brown
% % color_2 = [0.133, 0.545, 0.133]; % bright green
% 
% color_1 = [0, 0, 0]; % black
% color_2 = [1, 1, 1]; % white
% 
% cmap_custom = [linspace(color_1(1), color_2(1), Ncolors)', ...
%                linspace(color_1(2), color_2(2), Ncolors)', ...
%                linspace(color_1(3), color_2(3), Ncolors)'];
% 
% % Normalize RCS values for color mapping
% min_RCS = min(cell_RCS(:));
% max_RCS = max(cell_RCS(:));
% norm_RCS = (cell_RCS - min_RCS) / (max_RCS - min_RCS);
% 
% % Draw colored resolution cells
% for i = 1:n_ct_cells
%     for j = 1:n_at_cells
%         x_left = x_window(1) + (i-1)*res_ct;
%         x_right = x_left + res_ct;
%         y_bottom = y_window(1) + (j-1)*res_at;
%         y_top = y_bottom + res_at;
%         color_idx = round(norm_RCS(i,j) * (Ncolors-1)) + 1;
%         rect_color = cmap_custom(color_idx, :);
%         fill([x_left x_right x_right x_left], [y_bottom y_bottom y_top y_top], rect_color, 'EdgeColor', 'none');
%     end
% end
% 
% colormap(cmap_custom)
% cb = colorbar;
% cb.Label.String = 'Average RCS';  
% cb.Label.FontSize = 12;  
% cb.Label.FontWeight = 'bold';  
% cb.Label.Interpreter = 'latex';
% clim([min_RCS, max_RCS])
% axis equal
% axis([x_window y_window])
% set(gcf, 'Position', [80, 80, 600, 600])
% xlabel('Cross-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('Along-Track Distance [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% title('Simulated Ohia Radar Coverage', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
% set(gca, 'FontSize', 12, 'FontName', 'Arial');
% hold off

%% Other Plots

if export_std == 1
    exportgraphics(gcf, [folder,'\',fileName,'SAR','.pdf'], 'ContentType', 'vector')
end

if export_inf_chance == 1
    d_temp = 0:.1:30;
    figure()
    hold on
    plot(d_temp,inf_sprd(d_temp,1)*100,LineWidth=2)
    xlim([0 20])
    ylim([0 inf_sprd(0,1)*100])
    ylabel('Chance of Infection Spread [\%]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    xlabel('Distance from Host [m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'FontSize', 12, 'FontName', 'Arial');
    hold off
    exportgraphics(gcf, [folder,'\', fileName,'inf', '.pdf'],'ContentType','vector')
end

if export_dielecric == 1
    t_inf_temp = 0:.1:20; %days
    figure()
    hold on
    plot(t_inf_temp,leaf_vwc(t_inf_temp,dielectric_healthy,dielectric_infected,time_to_die),LineWidth=2)
    xlim([0 20])
    ylim([leaf_vwc(t_inf_temp(end),dielectric_healthy,dielectric_infected,time_to_die)-1 leaf_vwc(t_inf_temp(1),dielectric_healthy,dielectric_infected,time_to_die)+1])
    ylabel('Leaf Dielectric [F/m]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    xlabel('Time Infected [days]', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'FontSize', 12, 'FontName', 'Arial');
    hold off
    exportgraphics(gcf, [folder,'\', fileName, 'dielectric', '.pdf'],'ContentType','vector')
end
