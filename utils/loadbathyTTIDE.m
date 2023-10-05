function outputBathy = loadbathyTTIDE(lonlatlim, indstep)
% outputBathy = LOADBATHYTTIDE(lonlatlim, indstep)
%
%   inputs
%       - lonlatlim: 1x4 vector with longitude and latitude limits.
%       - indstep: an integer to skip data points.
%
%   outputs
%       - outputBathy: structure variable with bathymetry within lonlatlim.
%
%
%
% Olavo Badaro Marques.


%% Indicate bathymetry file path:

dirBathy = fullfile(paper_directory(), 'data', 'extradata');
        
% This file has the bathymetry of the Tasman slope/rise
% region (fields are lon/lat/depth):
bathyfile = 'TTIDE_bathymetry';

TTIDEbathy = load(fullfile(dirBathy, bathyfile));
TTIDEbathy = TTIDEbathy.ttide_bathymetry;


%%

if ~exist('lonlatlim', 'var') || isempty(lonlatlim)
    lonlim = [min(TTIDEbathy.lon) max(TTIDEbathy.lon)];
    latlim = [min(TTIDEbathy.lat) max(TTIDEbathy.lat)];
else
    lonlim = lonlatlim(1:2);
    latlim = lonlatlim(3:4);
end



%% Subset bathymetry:

subsetlon = TTIDEbathy.lon>=lonlim(1) & TTIDEbathy.lon<=lonlim(2);
subsetlat = TTIDEbathy.lat>=latlim(1) & TTIDEbathy.lat<=latlim(2);

TTIDEbathy.lon = TTIDEbathy.lon(subsetlon);
TTIDEbathy.lat = TTIDEbathy.lat(subsetlat);

TTIDEbathy.depth = TTIDEbathy.depth(subsetlat, subsetlon);
        

%%

if exist('indstep', 'var') 
    TTIDEbathy.lon = TTIDEbathy.lon(1:indstep:end);
    TTIDEbathy.lat = TTIDEbathy.lat(1:indstep:end);
    TTIDEbathy.depth = TTIDEbathy.depth(1:indstep:end, 1:indstep:end);
end


%%

outputBathy.lon = TTIDEbathy.lon;
outputBathy.lat = TTIDEbathy.lat;
outputBathy.depth = TTIDEbathy.depth;