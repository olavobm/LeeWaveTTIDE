function [xpts, ypts] = lonlat2kmgrid(lonlat0, lonpts, latpts)
% [xpts, ypts] = LONLAT2KMGRID(lonlat0, lonpts, latpts)
%
%   inputs
%       - lonlat0: 1x2 vector with the (longitude, latitude)
%                  coordinate reference (i.e. the 0 on the
%                  x/y grid).
%       - lonpts: array with longitudes of interest.
%       - latpts: array with latitudes, of the same size as lonpts.
%
%   outputs
%       - xpts: the x-coordinate of the inputs (relative to lonlat0).
%       - ypts: the y-coordinate.
%
%
%
% ONE OR MULITPLE POINTS???
%
% TO DO:
%   - Write comments.
%   - Include another input that allows for rotated coordinate system.
%
% See also: xy2lonlat.m
%
% Olavo Badaro Marques, 18/Jun/2018.



%%

[distpts, azymthpts] = distance(lonlat0(2), lonlat0(1), ...
                                 latpts, lonpts);

%
distpts_km = deg2km(distpts);


%%
%

%
xpts = distpts_km .* sind(azymthpts);

ypts = distpts_km .* cosd(azymthpts);