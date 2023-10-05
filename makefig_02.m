% Fig. 2. just the observation sites and the slope criticality.
%
%

clear
close all


%% Load bathymetry

%
lon_lims = [148.9, 149.08];
%
lat_lims = [-41.4, -41.25];

%
bathy = loadbathyTTIDE([lon_lims, lat_lims]);

%
bathy.depth = - double(bathy.depth);

% Make a grid
[bathy.longrid, bathy.latgrid] = meshgrid(bathy.lon, bathy.lat);


%% Bathymetry with criticality

%
load(fullfile(paper_directory(), 'data', 'extradata', 'TasmanBathy.mat'));

% Copy for the code block below
bathyAll = tasmanBathy;

%
lonlims_NTas = [148.645, 149.75];
latlims_NTas = [-41.4167, -41.2];

%
linlon_lims = (tasmanBathy.lon >= lonlims_NTas(1)) & ...
              (tasmanBathy.lon <= lonlims_NTas(2));
          
linlat_lims = (tasmanBathy.lat >= latlims_NTas(1)) & ...
              (tasmanBathy.lat <= latlims_NTas(2));
          
%
tasmanBathy.lon = tasmanBathy.lon(linlon_lims);
tasmanBathy.lat = tasmanBathy.lat(linlat_lims);
tasmanBathy.depth = tasmanBathy.depth(linlat_lims, linlon_lims);
tasmanBathy.slopeC = tasmanBathy.slopeC(linlat_lims, linlon_lims);
tasmanBathy.N2 = tasmanBathy.N2(linlat_lims, linlon_lims);

%
[tasmanBathy.longrid, ...
 tasmanBathy.latgrid] = meshgrid(tasmanBathy.lon, tasmanBathy.lat);


%% Same as above, but for the whole Tasman slope

%
lonlims_all = [147.7, 149.6];
latlims_all = [-43.6, -40.5];

%
linlon_all_lims = (bathyAll.lon >= lonlims_all(1)) & ...
                  (bathyAll.lon <= lonlims_all(2));
          
linlat_all_lims = (bathyAll.lat >= latlims_all(1)) & ...
                  (bathyAll.lat <= latlims_all(2));
          
%
bathyAll.lon = bathyAll.lon(linlon_all_lims);
bathyAll.lat = bathyAll.lat(linlat_all_lims);
bathyAll.depth = bathyAll.depth(linlat_all_lims, linlon_all_lims);
bathyAll.slopeC = bathyAll.slopeC(linlat_all_lims, linlon_all_lims);
bathyAll.N2 = bathyAll.N2(linlat_all_lims, linlon_all_lims);

%
[bathyAll.longrid, ...
 bathyAll.latgrid] = meshgrid(bathyAll.lon, bathyAll.lat);




%%

%
load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));

%
T1loc = mooringlocs.T1;
T2loc = mooringlocs.T2;


%% CTD/LADCP towyos

%
dir_shipdata = fullfile(paper_directory(), 'data', 'shipboard');

%
NT1data = load(fullfile(paper_directory(), 'data', 'shipboard', 'NT1', 'NT1_towyodata.mat'));

%
NT1data = NT1data.NT1data;


%% Fast-CTD towyo

%
FastCTDdata = load(fullfile(paper_directory(), 'data', 'shipboard', 'FastCTD', 'FastCTD_towyo.mat'));
FastCTDdata = FastCTDdata.FCTDtransect;


%% Shiptrack, to assign location to the LADCP data

%
ttide_shiptrack = load(fullfile(paper_directory(), 'data', 'metadata', 'TTIDE_shiptrack.mat'));
ttide_shiptrack = ttide_shiptrack.ttideshiptrack;

%
yday_station_lims = [45.95, 46.94];
%
linyday_lims = (ttide_shiptrack.yday >= yday_station_lims(1)) & ...
               (ttide_shiptrack.yday <= yday_station_lims(2));
indplt_lims = find(linyday_lims);

    
%%
% ---------------------------------------------------------------
% ------------ REFERENCE COORDINATE TO THE X/Y GRID -------------
% ---------------------------------------------------------------

%%

%
lonlat0 = [149.0019, -41.3348];


%% Reference bathymetry and observations to (x/y) grid

%
[bathy.xg, bathy.yg] = lonlat2kmgrid(lonlat0, bathy.longrid, bathy.longrid);


%% For the bathymetry with criticality

[tasmanBathy.xg, ...
 tasmanBathy.yg] = lonlat2kmgrid(lonlat0, tasmanBathy.longrid, ...
                                          tasmanBathy.latgrid);

                                      
%% T1 and T2 moorings

%
[T1loc.x, T1loc.y] = lonlat2kmgrid(lonlat0, T1loc.longitude, T1loc.latitude);
[T2loc.x, T2loc.y] = lonlat2kmgrid(lonlat0, T2loc.longitude, T2loc.latitude);


%% CTD/LADCP towyos
 
%
for i = 1:NT1data.nrepeats
	%
    [NT1data.CTDdata.transects(i).x, ...
     NT1data.CTDdata.transects(i).y] = ...
                    lonlat2kmgrid(lonlat0, ...
                                  NT1data.CTDdata.transects(i).lon, ...
                                  NT1data.CTDdata.transects(i).lat);
end


%% Fast-CTD         

%
for i = 1:FastCTDdata.Nrepeats
    
	%
    [FastCTDdata.timeseries(i).x, ...
     FastCTDdata.timeseries(i).y] = lonlat2kmgrid(lonlat0, ...
                                                  FastCTDdata.timeseries(i).lonstart, ...
                                                  FastCTDdata.timeseries(i).latstart);
end
 

%%

%
[ttide_shiptrack.x, ...
 ttide_shiptrack.y] = lonlat2kmgrid(lonlat0, ttide_shiptrack.lon, ...
                                             ttide_shiptrack.lat);
                                         

%%
% ---------------------------------------------------------------
% ------- GET BATHYMETRY AND CRITICALITY ALONG A TRANSECT -------
% ---------------------------------------------------------------

%% End points

% ------------------------------------------------------
% At the Fast-CTD transect line
transectref(1).lonwaypoints = [148.9584, 149.0327];
transectref(1).latwaypoints = -41.3265 .* [1, 1];
%
transectref(1).npts = 80;

% ------------------------------------------------------
% Mooring line line
transectref(2).lonwaypoints = [148.9584, 149.0327];
transectref(2).latwaypoints = -41.3349 .* [1, 1];
%
transectref(2).npts = 80;


% Extend further offshore
transectref(1).lonwaypoints(2) = 149.2;
transectref(2).lonwaypoints(2) = 149.2;


%%

%
for i = 1:length(transectref)
    
	%
    transectref(i).lonpoints = linspace(transectref(i).lonwaypoints(1), ...
                                        transectref(i).lonwaypoints(2), ...
                                        transectref(i).npts);

    transectref(i).latpoints = repmat(transectref(i).latwaypoints(1), 1, transectref(i).npts);

end


%% Create interpolants

bathy_interp = scatteredInterpolant(tasmanBathy.xg(:), ...
                                    tasmanBathy.yg(:), ...
                                    tasmanBathy.depth(:));   
                                
%
slopeC_interp = scatteredInterpolant(tasmanBathy.xg(:), ...
                                     tasmanBathy.yg(:), ...
                                     tasmanBathy.slopeC(:));   
          
                                 
%% Reference to x/y grid

%
for i = 1:length(transectref)
    
	%
    [transectref(i).x, ...
     transectref(i).y] = lonlat2kmgrid(lonlat0, transectref(i).lonpoints, ...
                                                transectref(i).latpoints);
                    
end

                                 
%%

%
for i = 1:length(transectref)
    
	%
    transectref(i).depth = bathy_interp(transectref(i).x, transectref(i).y);
    %
    transectref(i).slopeC = slopeC_interp(transectref(i).x, transectref(i).y);

end


%%
% -------------------------------------------------------------
% -------------------------------------------------------------
% ------------------------ MAKE FIGURE ------------------------
% -------------------------------------------------------------
% -------------------------------------------------------------


%
clrmap.gamma = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_velocity.mat'));
%
clrmap.gamma = clrmap.gamma.cmapredblue;


%%

% Get bathymetry and criticality averaged
% over some y-distance around the bump
indrowsavg = 120:1:149;
%
xavg_plt = nanmean(tasmanBathy.xg(indrowsavg, :), 1);
bathy_plt = nanmean(tasmanBathy.depth(indrowsavg, :), 1);
critic_plt = nanmean(tasmanBathy.slopeC(indrowsavg, :), 1);


%%

%
xaxs_sec_lim = [-4.7852, 4.5];
yaxs_sec_lim = [1500, 2250];

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.1527    0.1942    0.4232    0.4767];
%
haxs_2Dmap = axes('Position', [0.1, 0.3, 0.455, 0.455]);
haxs_section = axes('Position', [0.62, 0.3, 0.2, 0.455]);
haxs_crosssec = axes('Position', [0.775, 0.55, 0.175, 0.25]);
%
haxs_section.Position(1) = haxs_section.Position(1) - 0.025;
haxs_crosssec.Position(1) = haxs_crosssec.Position(1) - 0.025;
%
hold(haxs_2Dmap, 'on')
hold(haxs_section, 'on')
hold(haxs_crosssec, 'on')
    
% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------

    %
    pcolor(haxs_2Dmap, tasmanBathy.xg, tasmanBathy.yg, ...
                       log10(tasmanBathy.slopeC))
    %
    contour(haxs_2Dmap, tasmanBathy.xg, tasmanBathy.yg, ...
                        -tasmanBathy.depth, 1600:100:2400, ...
                        'k', 'LineWidth', 2)        
                    
%
shading(haxs_2Dmap, 'flat')
%
axis(haxs_2Dmap, 'equal')

% ----------------------------------------
%
statlbSZ = 86;

% T1
plot(haxs_2Dmap, T1loc.x, T1loc.y, '.k', 'MarkerSize', statlbSZ)

% T2
plot(haxs_2Dmap, T2loc.x, T2loc.y, '.k', 'MarkerSize', statlbSZ)

% NT1
for i = 1:NT1data.nrepeats
    %
    plot(haxs_2Dmap, NT1data.CTDdata.transects(i).x(:), ...
                     NT1data.CTDdata.transects(i).y(:), ...
                     'g', 'LineWidth', 4)
end

% Fast-CTD transect
for i = 1:FastCTDdata.Nrepeats
    %
    plot(haxs_2Dmap, FastCTDdata.timeseries(i).x(:), ...   
                     FastCTDdata.timeseries(i).y(:), ...
                     '-b', 'LineWidth', 2)
end

% Fast-CTD station
x_FCTD_station = mean(ttide_shiptrack.x(indplt_lims));
y_FCTD_station = mean(ttide_shiptrack.y(indplt_lims));
%
plot(haxs_2Dmap, x_FCTD_station, y_FCTD_station, '.b', 'MarkerSize', statlbSZ)
        
% ----------------------------------------
%
axes(haxs_2Dmap)

    %
    txtlblFS = 14;

    %
    ht_T1 = text(T1loc.x, T1loc.y, 'T1');
        ht_T1.FontSize = txtlblFS;
        ht_T1.Color = 'w';
        ht_T1.FontWeight = 'bold';
        ht_T1.HorizontalAlignment = 'center';

    %
    ht_T2 = text(T2loc.x, T2loc.y, 'T2');
        ht_T2.FontSize = txtlblFS;
        ht_T2.Color = 'w';
        ht_T2.FontWeight = 'bold';
        ht_T2.HorizontalAlignment = 'center';

    %
    ht_FCTD = text(x_FCTD_station, y_FCTD_station, 'F');
        ht_FCTD.FontSize = txtlblFS;
        ht_FCTD.Color = 'w';
        ht_FCTD.FontWeight = 'bold';
        ht_FCTD.HorizontalAlignment = 'center';

% NT1
fill([1.05, 1.05, 2.35, 2.35, 1.05], [1.25, 1.9, 1.9, 1.25, 1.25], 'k')
        
%
axes(haxs_2Dmap)

%
ht_NT1 = text(1.7112, 1.5833, 'NT1');
    ht_NT1.FontSize = txtlblFS;
    ht_NT1.Color = 'g';
    ht_NT1.FontWeight = 'bold';
    ht_NT1.HorizontalAlignment = 'center';

                
                
% ----------------------------------------
        
%
caxis(haxs_2Dmap, 0.8.*[-1, 1])

%
xaxs_lim = [-4.7852, 5.4743];
%         yaxs_lim = [-4.2193, 5.6952];
yaxs_lim = [-2.5, 5];
set(haxs_2Dmap, 'FontSize', 14, 'Box', 'on', ...
                'XLim', xaxs_sec_lim, 'YLim', yaxs_lim)

%
xlabel(haxs_2Dmap, '[km]', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel(haxs_2Dmap, '[km]', 'Interpreter', 'Latex', 'FontSize', 16)
%
title(haxs_2Dmap, 'Observation sites', 'Interpreter', 'Latex', 'FontSize', 20)


% ------------------------------------------------------------
% ------------------------------------------------------------
    

%
axes(haxs_section)

%
set(haxs_section, 'XLim', [-10, 10], 'YLim', [0, 4500])
%
mClrs_aux = matlabColors;
%
hf_subsec_aux = fill([xaxs_sec_lim(1), xaxs_sec_lim(1), xaxs_sec_lim(2), xaxs_sec_lim(2), xaxs_sec_lim(1)], ...
                     [2250, 1400, 1400, 2250, 2250], 'k');
    hf_subsec_aux.FaceColor = mClrs_aux(4, :);
    hf_subsec_aux.FaceAlpha = 0.5;
    hf_subsec_aux.LineStyle = 'none';


%
hfBathy = fill2Dbathy(xavg_plt, -bathy_plt);
    hfBathy.FaceColor = 0.5.*[1, 1, 1];

%
%         scatter(xavg_plt, -bathy_plt, 10, log10(critic_plt), 'filled')
scatter(xavg_plt, -bathy_plt, 10, log10(critic_plt))
            
      

% --------------------------------

%
transectSlopes(1).xlim = [-28.45, -11];
%     transectSlopes(2).xlim = [-7.42, 10.7];
transectSlopes(2).xlim = [1, 10];
transectSlopes(3).xlim = [11.0860, 23.03];

%
for i = 1:length(transectSlopes)
    %
    linxlims_aux = (xavg_plt >= transectSlopes(i).xlim(1)) & ...
                   (xavg_plt <= transectSlopes(i).xlim(2));

    %
    transectSlopes(i).zlim = -[bathy_plt(find(linxlims_aux, 1, 'first')), ...
                               bathy_plt(find(linxlims_aux, 1, 'last'))];

    %
    transectSlopes(i).criticavg = nanmean(critic_plt(linxlims_aux));
end

%
plot(haxs_section, transectSlopes(1).xlim, ...
                   transectSlopes(1).zlim+280, '-.k', 'LineWidth', 2)
%
plot(haxs_section, transectSlopes(2).xlim, ...
                   transectSlopes(2).zlim+300, '-.k', 'LineWidth', 2)
%
plot(haxs_section, transectSlopes(3).xlim, ...
                   transectSlopes(3).zlim+130, '-.k', 'LineWidth', 2)

%
htxt_1_aux = text(-22.5, 1550, num2str(transectSlopes(1).criticavg, '%.1f'));
htxt_2_aux = text(-13, 2800, ['$\gamma$ = ' num2str(transectSlopes(2).criticavg, '%.1f')]);
htxt_3_aux = text(13, 3850, num2str(transectSlopes(3).criticavg, '%.1f'));
      
%
txt_FS = 18;
htxt_1_aux.FontSize = txt_FS;
htxt_2_aux.FontSize = txt_FS;
htxt_3_aux.FontSize = txt_FS;

%
htxt_1_aux.HorizontalAlignment = 'center';
htxt_2_aux.HorizontalAlignment = 'center';
htxt_3_aux.HorizontalAlignment = 'center';

%
htxt_1_aux.Interpreter = 'Latex';
htxt_2_aux.Interpreter = 'Latex';
htxt_3_aux.Interpreter = 'Latex';

    
% --------------------------------

%
hcb_2Dmap = colorbar(haxs_section);

    %
    hcb_2Dmap.FontSize = 12;
    hcb_2Dmap.Location = 'eastoutside';
    hcb_2Dmap.Position = [0.8, 0.3, 0.0125, 0.2];
    hcb_2Dmap.Label.Interpreter = 'Latex';
    hcb_2Dmap.Label.String = '$\log_{10} \gamma$';
    hcb_2Dmap.Label.FontSize = 18;
    hcb_2Dmap.Label.Rotation = 0;
    hcb_2Dmap.Label.VerticalAlignment = 'middle';
    hcb_2Dmap.Label.Position = [6.75, 0, 0];
    hcb_2Dmap.Ticks = -0.6:0.3:0.6;


%
set(haxs_section, 'FontSize', 12, 'Box', 'on', ...
                             'YDir', 'reverse', ...
                             'CLim', [-0.8, 0.8], ...
                             'XLim', [-29, 50], 'YLim', [0, 4200])
set(haxs_section, 'YTick', 0 : 1000 : 4000)
set(haxs_section, 'Color', 0.9.*[1, 1, 1])

%
xlabel(haxs_section, '[km]', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel(haxs_section, 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 16)

    
% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------
%
indtrsctplt = 2;

%
axes(haxs_crosssec)

    %
    set(haxs_crosssec, 'XLim', xaxs_sec_lim, ...
                                  'YLim', yaxs_sec_lim)
    overlayline(haxs_crosssec, 'v', [T1loc.x, T2loc.x], '-k', 'LineWidth', 2)
    overlayline(haxs_crosssec, 'v', mean(ttide_shiptrack.x(indplt_lims)), '-b', 'LineWidth', 2)


%
hfBathy = fill2Dbathy(xavg_plt, -bathy_plt);
    hfBathy.FaceColor = 0.5.*[1, 1, 1];
    
    %
    scatter(xavg_plt, -bathy_plt, 70, log10(critic_plt), 'filled');

%
caxis(haxs_crosssec, caxis(haxs_2Dmap))

% -----------------------------------
% Add values of criticality
criticchunks(1).xlims = [-0.8748, -0.0899];
criticchunks(1).linlims = (xavg_plt >= criticchunks(1).xlims(1)) & ...
                          (xavg_plt <= criticchunks(1).xlims(2));
criticchunks(1).meancritic = mean(critic_plt(criticchunks(1).linlims));
criticchunks(1).zlims = -bathy_plt(criticchunks(1).linlims);
criticchunks(1).zlims = criticchunks(1).zlims([1, end]);
%
criticchunks(2).xlims = [1.2556, 2.8814];
criticchunks(2).linlims = (xavg_plt >= criticchunks(2).xlims(1)) & ...
                          (xavg_plt <= criticchunks(2).xlims(2));
criticchunks(2).meancritic = mean(critic_plt(criticchunks(2).linlims));
criticchunks(2).zlims = -bathy_plt(criticchunks(2).linlims);
criticchunks(2).zlims = criticchunks(2).zlims([1, end]);
%
criticchunks(3).xlims = [0.1063, 1.1715];
criticchunks(3).linlims = (xavg_plt >= criticchunks(3).xlims(1)) & ...
                          (xavg_plt <= criticchunks(3).xlims(2));
criticchunks(3).meancritic = mean(critic_plt(criticchunks(3).linlims));
criticchunks(3).zlims = -bathy_plt(criticchunks(3).linlims);
criticchunks(3).zlims = criticchunks(3).zlims([1, end]);

%
indchunkplt = [1, 2];
for i = indchunkplt
    plot(haxs_crosssec, criticchunks(i).xlims, ...
                        criticchunks(i).zlims + 80, '-.k', 'LineWidth', 2)
end
        
%
axes(haxs_crosssec)
    %
    xtxt_1 = -3.8;
    ytxt_1 = 2100;
    htxt_1 = text(xtxt_1, ytxt_1, ...
                  ['$\gamma$=' num2str(criticchunks(1).meancritic, '%.1f')]);
    htxt_1.Interpreter = 'Latex';
    htxt_1.FontSize = 20;
    %
    xtxt_2 = 1.25;
    ytxt_2 = 2200;
    htxt_2 = text(xtxt_2, ytxt_2, ...
                  [num2str(criticchunks(2).meancritic, '%.1f')]);
% %                           ['$\gamma$=' num2str(criticchunks(2).meancritic, '%.1f')]);

    htxt_2.Interpreter = 'Latex';
    htxt_2.FontSize = 20;

    %
    lblmoorFS = 16;
    %
    htxt_T1moor = text(-3, 1600, 'T1');
    htxt_T2moor = text(1.6, 1600, 'T2');
    htxt_Fstatn = text(-0.6, 1580, 'F');
    %
    htxt_T1moor.FontSize = lblmoorFS;
    htxt_T2moor.FontSize = lblmoorFS;
    htxt_Fstatn.FontSize = lblmoorFS;
    %
    htxt_Fstatn.Color = [0, 0, 1];
            
% -----------------------------------
%
colormap(haxs_2Dmap, clrmap.gamma);
colormap(haxs_section, clrmap.gamma);
colormap(haxs_crosssec, clrmap.gamma);

%
set(haxs_crosssec, 'FontSize', 12, 'Box', 'on', ...
                   'XGrid', 'on', 'YGrid', 'on', ...
                   'XLim', xaxs_sec_lim, ...
                   'YLim', yaxs_sec_lim, 'YDir', 'reverse')
set(haxs_crosssec, 'YTick', 1600:200:2200, ...
                   'YTickLabel', {'', '1800', '2000', '2200'})

%
set(haxs_crosssec, 'Color', 0.9 .* [1, 1, 1])

    
% ------------------------------------------------------------

%
axes(haxs_2Dmap)

    %
    lettersFS = 16;
    yheight = 5.5;

    %
    htxt_a = text(-4.8, yheight, 'a)');
        htxt_a.FontSize = lettersFS;
    %
    htxt_b = text(7.2, yheight, 'b)');
        htxt_b.FontSize = lettersFS;

    %
    htxt_title = text(7.9, yheight, 'Criticality');
        htxt_title.FontSize = 20;
        htxt_title.Interpreter = 'Latex';


%
set(haxs_2Dmap, 'CLim', 0.65.*[-1, 1])
    
    
%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure02.pdf'), 'Resolution', 300)

