clear;close all;clc;
% Define your parameters
f_param = 1.6;          % Control overall smoothness of noise
initialH = 3;          % Base height
initialP = 20;          % Initial overall height variation
min_X = 0; max_X = 100;
min_Y = 0; max_Y = 100;
num_Iter = 7;           % More iterations for denser perturbations
valley_exp = 0.5;       % Key for wider peaks, shallower valleys (0 to 1)

% Generate the terrain
[x, y, terrain] = helperRandomTerrainGenerator(f_param, initialH, ...
    initialP, min_X, max_X, min_Y, max_Y, num_Iter, valley_exp);

% Plot the terrain (assuming helperPlotSimulatedTerrain is in your path)
helperPlotSimulatedTerrain(x, y, terrain);

function [x,y,terrain] = helperRandomTerrainGenerator(f,initialHeight, ...
    initialPerturb,minX,maxX,minY,maxY,numIter, valley_exponent)
%helperRandomTerrainGenerator Generate forest canopy-like terrain
% [x,y,terrain] = helperRandomTerrainGenerator(f,initialHeight, ...
%    initialPerturb,minX,maxX,minY,maxY,numIter, valley_exponent)
%
% Inputs:
%   - f                   = A roughness parameter. Higher values lead to
%                           smoother terrain with less impact from finer
%                           perturbations.
%   - initialHeight       = Sets the initial height of the lattice.
%   - initialPerturb      = Initial perturbation amount.
%   - minX,maxX,minY,maxY = Initial points defining the extent.
%   - numIter             = Number of iterations, affecting detail density.
%   - valley_exponent     = A value > 1 for pushing low values up (shallower valleys).
%                           Typical range: 1.5 to 4.
%
% Output:
%    - x                  = X-dimension mesh grid
%    - y                  = Y-dimension mesh grid
%    - terrain            = Two-dimensional array representing terrain/canopy height.
%
% Generate random terrain with wider, rounded peaks and shallower, rarer valleys.

dX = (maxX-minX)/2;
dY = (maxY-minY)/2;
[x,y] = meshgrid(minX:dX:maxX,minY:dY:maxY);
terrain = ones(3,3)*initialHeight;
perturb = initialPerturb;

% Define a Gaussian smoothing filter (good for rounding)
% Adjust kernelSize and sigma for desired overall smoothness of features.
kernelSize = [7 7]; % Larger kernel for more rounding
sigma = 1.8;        % Larger sigma for more blurring/smoothing
smoothingKernel = fspecial('gaussian', kernelSize, sigma);

for ii = 2:numIter
    perturb = perturb/f;
    
    oldX = x;
    oldY = y;
    dX = (maxX-minX)/2^ii;
    dY = (maxY-minY)/2^ii;
    [x,y] = meshgrid(minX:dX:maxX,minY:dY:maxY);
    terrain = griddata(oldX,oldY,terrain,x,y);
    
    terrain = terrain + perturb * random('norm',0,1,size(x,1),size(x,2));
    
    % Apply smoothing to round off peaks and valleys
    terrain = imfilter(terrain, smoothingKernel, 'replicate');
    
    % --- NEW: Height Remapping to control peaks and valleys ---
    % Normalize terrain to a [0, 1] range temporarily for remapping
    min_val = min(terrain(:));
    max_val = max(terrain(:));
    
    if (max_val - min_val) > 1e-6 % Avoid division by zero for flat terrain
        normalized_terrain = (terrain - min_val) / (max_val - min_val);
        
        % Apply power curve: x^exponent.
        % If exponent > 1, it pushes lower values UP and spreads higher values.
        % This makes valleys shallower and peaks wider/more prominent.
        remapped_terrain = normalized_terrain .^ valley_exponent;
        
        % Scale back to original height range (or a desired range)
        terrain = min_val + remapped_terrain * (max_val - min_val);
    end
    % --- END NEW ---

    % Ensure no negative heights (after remapping)
    terrain(terrain < 0) = 0;
end
end

function helperPlotSimulatedTerrain(xvec,yvec,A)
% Plot simulated terrain

figure()
hS = surf(xvec,yvec,A);
hS.EdgeColor = 'none';
hC = colorbar;
hC.Label.String = 'Elevation (m)';
landmap = landColormap(64);
colormap(landmap); 
xlabel('X (m)')
ylabel('Y (m)')
axis equal;
title('Simulated Terrain')
view([78 78])
drawnow
pause(0.25)
end