function [pt_tree,N_trees,x_window,y_window] = tree_gen(n_ct_cells,n_at_cells,res_ct,res_at,density_ohia,percent_infected,ohia_diameters)
% Generates trees and grid for Ohia Sim.
% Inputs: 
%       n_at/ct_cells = #along track and cross treck cells
%       density_ohia = tree density in #/area
%       percent_infected = $ of Ohia trees to have ROD
%       ohia_diameters = [min_avg diameter [m], [max_avg diameter [m]]
% Output:
%       pt_tree = [x_location, y_location, is_infected (bool), diameter]


    
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
    pt_tree(:,3) = rand(N_trees,1) < percent_infected; % 1 if infected, 0 if healthy
    pt_tree(:,4) = (ohia_diameters(1) + (ohia_diameters(2)-ohia_diameters(1)).*rand(N_trees,1))/2;
    pt_tree(:,5) = 0; %time infected
    
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
end