function cell_RCS = get_cell_RCS(x_window, y_window, n_ct_cells, n_at_cells, res_ct, res_at, sample_distance, RCS_base, pt_tree, dielectric_healthy, dielectric_infected, dielectric_std, wavelength,time_to_die)
%CALCULATE_RCS Calculates the Radar Cross Section (RCS) for each resolution cell.
%
%   This function calculates the RCS for each cell in a grid, considering the
%   presence of trees (healthy and infected) within each cell.
%
%   Input Arguments:
%   x_window, y_window:  Vectors defining the spatial extent of the area of interest [m].
%   n_ct_cells, n_at_cells: Number of cells in the cross-track and along-track directions.
%   res_ct, res_at:    Cross-track and along-track resolution [m].
%   sample_distance:   Distance between samples within each resolution cell [m].
%   RCS_base:          Base RCS value for the ground [dB].
%   pt_tree:            Matrix containing tree information [x, y, is_infected, radius].
%   dielectric_healthy, dielectric_infected: Dielectric constants for healthy and infected trees.
%   dielectric_std:    Standard deviation of the dielectric constant.
%   wavelength:        Wavelength of the radar signal [m].
%
%   Output Arguments:
%   cell_RCS:          A matrix containing the average RCS value for each
%                      resolution cell [dB].
%
%   Detailed Explanation:
%   1.  Creates a grid of resolution cells based on the input parameters.
%   2.  Generates sample points within each resolution cell.
%   3.  Iterates through each resolution cell and calculates the RCS by:
%       a.  Initializing the RCS of the cell to the base ground RCS.
%       b.  Checking each sample point within the cell against all trees.
%       c.  If a sample point falls within a tree's radius, the RCS of that
%           sample is updated based on the tree's dielectric properties
%           (healthy or infected).  The function assumes the radar signal
%           cannot penetrate through multiple trees or the ground, adding
%           the RCS only once.
%       d.  If a sample point is not within any tree, its RCS is averaged
%           with the rest of the cell.
%       e.  The average RCS of all samples within the cell is assigned to
%           the cell's RCS value.
%   4. Returns the cell_RCS matrix.
%
%   Assumptions:
%   -   Trees are modeled as circles.
%   -   RCS values from trees are added linearly (not truly accurate, but a simplification).
%   -   Each resolution cell's RCS is the average of its samples.
%

    % Cell Creation
    cell_RCS = RCS_base*ones(n_ct_cells,n_at_cells); % Preallocate RCS grid
    N_trees = size(pt_tree, 1);


    % Create samples inside each cell
    samples_ct_cell = 0:sample_distance:res_ct;
    samples_at_cell = 0:sample_distance:res_at;
    [X, Y] = meshgrid(samples_ct_cell, samples_at_cell);
    xy_pairs = [X(:)+x_window(1), Y(:)+y_window(1)];
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
                for m = N_trees:-1:1
                    %cycle backward here because the tree "on top" is the last
                    %plotted tree, so have to reverse
                    x_tree = pt_tree(m,1);
                    y_tree = pt_tree(m,2);
                    is_infected = pt_tree(m,3);
                    r_tree = pt_tree(m,4);
                    d = sqrt((x_sample-x_tree)^2 + (y_sample-y_tree)^2);
                    if d < r_tree
                        %when adding the RCS value, the radar cannot see
                        %through multiple overlapping trees and the ground, so
                        %just add once
                        if RCS == RCS_base
                            if is_infected == 1
                                t_inf = pt_tree(m,5);
                                % dielectric = normrnd(dielectric_infected, dielectric_std);
                                dielectric = leaf_vwc(t_inf,dielectric_healthy,dielectric_infected, time_to_die);
                            else
                                dielectric = normrnd(dielectric_healthy, dielectric_std);
                            end
                            n = sqrt(dielectric);
                            k_squared = abs((n^2-1)/(n^2+1));
                            sigma_0 = pi^5*k_squared/wavelength^4;
                            RCS = sigma_0;    % Linear
                            RCS = 10*log10(RCS); %dB
                        end
                    end
                end
                samples_cell(l,3) = RCS;
                if samples_cell(l,3) == RCS_base % If its an empty spot
                    samples_cell(l,3) = mean(cell_RCS(:));%just average with rest
                end
            end
            cell_RCS(i,j) = mean(samples_cell(:,3)); % Average RCS in cell
            samples_cell(:,2) = samples_cell(:,2) + res_at; % Shift Y
        end
        samples_cell(:,2) = samples_cell(:,2) - n_at_cells*res_at;
        samples_cell(:,1) = samples_cell(:,1) + res_ct; % Shift X
    end
end
