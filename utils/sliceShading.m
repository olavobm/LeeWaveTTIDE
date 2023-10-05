function sliceShading(ts, clrs, alphatransp, tjoin)
% SLICESHADING(ts, clrs, alphatransp, tjoin)
%
%   inputs:
%       - ts: Nx2, 
%       - clrs: Nx3, N colors.
%       - alphatransp: from 0 to 1 (0 = fully transparent). Default is 0.5.
%
% Function SLICESHADING.m adds shading to the current time series plot.
%
% Include tjoin (the time limit which would join two time intervals
% in one shading if the intervals are separated at most by tjoin).
%
% Suggestions for the future:
%   - Later on, implements color choosing by the piece
%     of code inside stackedSahding.m
%   - Would also be nice to include labels for each shaded area.
%   - I can eventually extend this for some fancy 3D application.
%   - control alpha/transparency from input.
%
% Olavo Badaro Marques, 08/Dec/2016.


%%

if ~exist('alphatransp', 'var')
	alphatransp = 0.5;
end

if (alphatransp>0) && (alphatransp<0.05)
    error('I have updated sliceShading.m. Change the function call.')
end


%%

if isvector(clrs)
	clrs = repmat(clrs, size(ts, 1), 1);
end


%% Create vector of the y coordinates of the 

ylimsplt = ylim;
ylimsplt = [ylimsplt(1); ylimsplt(1); ylimsplt(2); ylimsplt(2)];


%% Create a 4xN array, where each column are the x
% coordinate of the vertices to be shaded:

tsaux = ts';
tsplt = [tsaux; flipud(tsaux)];


%% Loop through the columns of tsplt and shaded each slice:

hold on

%
for i = 1:size(ts, 1)
    
    hfaux = fill(tsplt(:, i), ylimsplt, clrs(i, :));
    
        hfaux.EdgeColor = 'none';  % because there is no transparency for the edge
        hfaux.FaceAlpha = alphatransp;
        
    %     set(hdpaux, 'LineStyle', 'none')
    
end


