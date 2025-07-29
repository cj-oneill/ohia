function [x, y, terrain] = helperForestCanopyGenerator(numIter, xmin, xmax, ymin, ymax, roughness, initialHeight, initialPerturb, peakSharpness, smoothing)
%helperForestCanopyGenerator Generates a random terrain surface resembling a forest canopy.
%
%   This function uses an iterative fractal algorithm (midpoint displacement)
%   to build a terrain. It includes parameters for controlling the overall
%   shape and texture to simulate the bumpy, overlapping surface of treetops.
%
% SYNTAX:
%   [x, y, terrain] = helperForestCanopyGenerator(numIter, xmin, xmax, ymin, ymax, roughness, initialHeight, initialPerturb, peakSharpness, smoothing)
%
% INPUTS:
%   - numIter:         Number of iterations. This controls the detail level and
%                      grid density. The final grid will be (2^n+1 x 2^n+1).
%                      A value between 5 and 9 is typically good.
%
%   - xmin, xmax:      The minimum and maximum bounds for the X-axis.
%
%   - ymin, ymax:      The minimum and maximum bounds for the Y-axis.
%
%   - roughness:       A roughness parameter. Lower values (e.g., 1.5-2.0)
%                      result in a rougher, more "bumpy" terrain. Higher
%                      values create a smoother surface.
%
%   - initialHeight:   The base height of the canopy layer before perturbations
%                      are added.
%
%   - initialPerturb:  The initial magnitude of random height variations. This
%                      sets the overall height difference between the tops of
%                      the trees and the gaps between them.
%
%   - peakSharpness:   An exponent for remapping heights. This is key for the
%                      canopy effect.
%                      - Values < 1.0 (e.g., 0.6-0.8) create shallower valleys
%                        and more prominent, rounded peaks, like treetops.
%                      - A value of 1.0 has no effect.
%                      - Values > 1.0 create sharper valleys and flatter peaks.
%
%   - smoothing:       The standard deviation (sigma) of the Gaussian smoothing
%                      kernel applied at each iteration. This rounds the sharp
%                      edges to create a more organic, "treetop" look. A value
%                      of 1 to 2 is a good starting point.
%
% OUTPUTS:
%   - x:               X-coordinates of the terrain mesh.
%   - y:               Y-coordinates of the terrain mesh.
%   - terrain:         A matrix where each value is the height of the canopy.
%
% EXAMPLE USAGE:
%   % --- Parameters for a dense, relatively uniform canopy ---
%   numIter = 8;         % Corresponds to 2^8+1 = 257 grid size
%   xmin = 1000; xmax = 1256; % Define X bounds
%   ymin = 1000; ymax = 1256; % Define Y bounds
%   roughness = 1.8;     % Lower value for more bumpiness
%   initialHeight = 50;  % Base height of the forest
%   initialPerturb = 20; % Tree heights vary by about +/- 20
%   peakSharpness = 0.7; % Make valleys shallower to represent dense canopy
%   smoothing = 1.5;     % Smooth the tops of the trees
%
%   % --- Generate the terrain ---
%   [x, y, terrain] = helperForestCanopyGenerator(numIter, xmin, xmax, ymin, ymax, ...
%       roughness, initialHeight, initialPerturb, peakSharpness, smoothing);
%
%   % --- Visualize the result ---
%   figure;
%   s = surf(x, y, terrain);
%   s.EdgeColor = 'none';
%   s.FaceLighting = 'gouraud';
%   axis equal;
%   title('Generated Forest Canopy');
%   xlabel('X'); ylabel('Y'); zlabel('Height');
%   colormap(summer); % Green colormap for a forest look
%   camlight head;    % Add a light source
%   material dull;

% --- 1. Input Validation and Initialization ---

if nargin < 5
    error('Not enough input arguments. Please provide numIter, xmin, xmax, ymin, and ymax.');
end

% Create an initial small grid (3x3) using the specified bounds.
[x, y] = meshgrid(linspace(xmin, xmax, 3), linspace(ymin, ymax, 3));
terrain = ones(3, 3) * initialHeight;

% The perturbation amount will decrease with each iteration.
perturb = initialPerturb;

% Define the smoothing kernel once.
smoothingKernel = fspecial('gaussian', [7 7], smoothing);

% --- 2. Iterative Terrain Generation ---
for ii = 2:numIter
    
    % Reduce perturbation amount for finer details.
    perturb = perturb / roughness;
    
    % Store the old grid and create the new, finer grid.
    oldX = x;
    oldY = y;
    oldTerrain = terrain;
    
    currentGridPoints = 2^ii + 1;
    [x, y] = meshgrid(linspace(xmin, xmax, currentGridPoints), linspace(ymin, ymax, currentGridPoints));
    
    % Interpolate the old, coarse terrain onto the new, finer grid.
    % 'linear' is fast and works well here.
    terrain = interp2(oldX, oldY, oldTerrain, x, y, 'linear');
    
    % Add random noise to the new grid.
    % The noise magnitude is controlled by 'perturb'.
    noise = random('norm', 0, 1, currentGridPoints, currentGridPoints);
    terrain = terrain + perturb * noise;
    
    % Apply smoothing to round the features, making them look more like treetops.
    terrain = imfilter(terrain, smoothingKernel, 'replicate');
end

% --- 3. Height Remapping for Canopy Effect ---
% This step is crucial for controlling the overall shape of the landscape.
min_val = min(terrain(:));
max_val = max(terrain(:));

% Avoid division by zero if terrain is completely flat.
if (max_val - min_val) > 1e-6
    % Normalize terrain to a [0, 1] range.
    normalized_terrain = (terrain - min_val) / (max_val - min_val);
    
    % Apply the power curve using the peakSharpness exponent.
    % An exponent < 1 pushes mid-range values UP towards the max value.
    % This makes the valleys shallower and the peaks more prominent, which
    % is characteristic of a dense forest canopy where gaps don't go to the ground.
    remapped_terrain = normalized_terrain .^ peakSharpness;
    
    % Scale the remapped terrain back to its original height range.
    terrain = min_val + remapped_terrain * (max_val - min_val);
end

% --- 4. Final Touches ---
% Ensure the canopy doesn't go below a certain height (e.g., the base height).
terrain(terrain < initialHeight-1) = initialHeight;

end
