function [B, circleIndexCounter] = circleGen(x_grid,y_grid,A,radiusMean,radiusStdDev)
    % Get rows and columns from the generated A
    [rows, cols] = size(A);

    % --- Initialize Matrices ---
    coveredCells = false(rows, cols);
    B = zeros(rows, cols);
    drawnCirclesData = {};
    circleIndexCounter = 0;

    % --- Nested Function: boxMullerTransform (Normal Distribution) ---
    % Generates a random number from a normal distribution.
    % mean: The mean of the distribution.
    % stdDev: The standard deviation of the distribution.
    function r = boxMullerTransform(mean, stdDev)
        u = 0; v = 0;
        % Ensure u and v are not zero to avoid log(0)
        while u == 0, u = rand(); end
        while v == 0, v = rand(); end
        z = sqrt(-2.0 * log(u)) * cos(2.0 * pi * v);
        r = z * stdDev + mean;
        % Ensure radius is positive and has a minimum value
        r = max(0.01, r); % Minimum radius of 0.01 meters
    end

    % --- Main Logic to Draw Circles ---
    fprintf('Starting circle generation...\n');
    totalCells = rows * cols;
    fprintf('Circles to draw: %i\n', totalCells)
    coveredCount = 0;
    iteration = 0;
    maxIterations = totalCells * 2; % Safety break to prevent infinite loops

    while coveredCount < totalCells && iteration < maxIterations
        % --- 1. Find the highest uncovered point ---
        maxVal = -inf; % Initialize with negative infinity for finding max
        highestR = -1; % Row index of the highest uncovered point
        highestC = -1; % Column index of the highest uncovered point

        for r = 1:rows
            for c = 1:cols
                % MATLAB uses 1-based indexing
                % Check if the cell is not covered and its height is greater than current max
                if ~coveredCells(r, c) && A(r, c) > maxVal
                    maxVal = A(r, c);
                    highestR = r;
                    highestC = c;
                end
            end
        end

        if highestR == -1
            % No more uncovered points found, break the loop
            fprintf('No more uncovered points found. Breaking loop.\n');
            break;
        end

        % Increment circle counter
        circleIndexCounter = circleIndexCounter + 1;

        % --- 2. Generate radius using normal distribution (this is in meters) ---
        currentRadiusMeters = boxMullerTransform(radiusMean, radiusStdDev);

        % Store circle data (optional, for analysis)
        drawnCirclesData{end+1} = struct('row', highestR, 'col', highestC, ...
                                         'radius', currentRadiusMeters, 'height', maxVal, ...
                                         'index', circleIndexCounter);

        % --- 3. Mark cells covered by the new circle and assign circle value to B ---
        % Iterate through all cells to check if they are within the current circle's radius
        for r_idx = 1:rows
            for c_idx = 1:cols
                % Calculate the physical (meter) coordinates of the current cell and the circle center
                % Use the x_grid and y_grid matrices directly for meter coordinates
                current_x_meter = x_grid(r_idx, c_idx);
                current_y_meter = y_grid(r_idx, c_idx);
                center_x_meter = x_grid(highestR, highestC);
                center_y_meter = y_grid(highestR, highestC);

                % Calculate distance in meters from current cell to circle center
                dist_meters = sqrt((current_x_meter - center_x_meter)^2 + (current_y_meter - center_y_meter)^2);

                % If the cell is within the circle's radius (in meters)
                if dist_meters <= currentRadiusMeters
                    % If the cell is not yet covered, mark it and assign the current circle's index
                    if ~coveredCells(r_idx, c_idx)
                        coveredCells(r_idx, c_idx) = true;
                        B(r_idx, c_idx) = circleIndexCounter; % Assign the circle's index
                    end
                end
            end
        end

        % Recalculate coveredCount
        coveredCount = sum(coveredCells(:));
        iteration = iteration + 1;

        % fprintf('Iteration %d: Covered %d/%d cells. Circles drawn: %d\n', ...
        %         iteration, coveredCount, totalCells, circleIndexCounter);
    end

    fprintf('Circle generation complete.\n');
    fprintf('Total circles drawn: %d\n', circleIndexCounter);

    % 
    % % --- Visualize Matrix B (simple heatmap) ---
    % figure;
    % % Use x_grid(1,:) and y_grid(:,1) for axes, as they represent the meter ranges
    % imagesc(x_grid(1,:), y_grid(:,1), B);
    % colorbar;   % Show color scale
    % colormap(jet); % Using 'jet' for a more traditional rainbow
    % axis xy; % Ensure y-axis is oriented correctly for imagesc
    % axis equal tight; % Keep aspect ratio and remove extra space
    % title('Matrix B: Lowest Circle Index for Each Cell');
    % xlabel('X (meters)');
    % ylabel('Y (meters)');

    % --- Plot B on the 3D surface of A ---
    figure; % Create a new figure for the 3D plot
    hold on
    num_unique_circles = max(B(:));
    if num_unique_circles == 0
        num_unique_circles = 1; % Avoid division by zero if no circles are drawn
    end
    num_colors_per_cycle = 64; % Standard number of colors in 'jet'
    num_cycles = 1;%ceil(num_unique_circles / num_colors_per_cycle) + 1; % Added +1 for better repetition

    % Create a repeated colormap
    repeated_colormap = repmat(jet(num_colors_per_cycle), num_cycles, 1);


    surf(x_grid, y_grid, A, B,EdgeColor="none"); % Plot A as a surface, with colors determined by B, using x and y for axes
    zlim([0 max(max(A))+1])
    colorbar;   % Show color scale for B values
    colormap(repeated_colormap); % Apply the custom repeated colormap
    title('Individual Tree Map Overlay');
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Height (meters)');
    view(3); % Set the default 3D view
    axis tight; % Adjust axis limits to tightly fit the data
    x_range = max(max(x_grid))-min(min(x_grid));
    y_range = max(max(y_grid))-min(min(y_grid));
    z_range = max(max(A)) - 0;
    zlim([0 max(max(A))])
    yx_ratio = y_range/x_range;
    zx_ratio = 2;
    zlim([0 max(max(A))])
    daspect([x_range y_range/yx_ratio z_range*zx_ratio]);    hold off
end
