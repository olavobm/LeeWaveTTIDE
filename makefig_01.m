%% Fig. 1. map of internal tide and dissipation in the model.

clear
close all



%% Load data

% Mooring locations
moorlocs = load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));
moorlocs = moorlocs.mooringlocs;

% Dissipation
modelD = load(fullfile(paper_directory(), 'data', 'extradata', 'data_modeldissipation.mat'));

% Altimetry
altdata = load(fullfile(paper_directory(), 'data', 'extradata', 'data_altimetry.mat'));


%% Convert x/y grid to lon/lat

% Convert from m to km (because JK_xy2lonlat.m need x,y in km)
modelD.x = modelD.x./1000;
modelD.y = modelD.y./1000;
            
%
modelD = JK_xy2lonlat(modelD);
    

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


%%

%
mooringlocs = load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));
mooringlocs = mooringlocs.mooringlocs;


%%
% -------------------------------------------------------
% --------------------- MAKE FIGURE ---------------------
% -------------------------------------------------------


%% Load colormaps

%
clrmap.dissipation = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_epsilon.mat'));
clrmap.ssh = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_velocity.mat'));

%
clrmap.dissipation = clrmap.dissipation.cmaporange;
clrmap.ssh = clrmap.ssh.cmapredblue;


%%

%
lon_lims_plt = [147.8, 149.5];
lat_lims_plt = [-43.5, -40.8];

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.2, 0.2, 0.7, 0.6];

%
haxs_alt = axes('Position', [0.4, 0.2, 0.5, 0.5]);
haxs_dsp = axes('Position', [0.05, 0.2, 0.5, 0.5]);
%
hold(haxs_alt, 'on');
hold(haxs_dsp, 'on');
    
% ---------------------------------------------
% ---------------------------------------------
%
axes(haxs_alt)

    %
    lonlims_sub = [146.3, 166.7];
    latlims_sub = [-49.8, -41];
    %
    linlonlims_aux = (altdata.M2_Tasman.SUM.lon >= lonlims_sub(1)) & ...
                     (altdata.M2_Tasman.SUM.lon <= lonlims_sub(2));
    linlatlims_aux = (altdata.M2_Tasman.SUM.lat >= latlims_sub(1)) & ...
                     (altdata.M2_Tasman.SUM.lat <= latlims_sub(2));

    %
    m_proj('lambert', 'long', [146, 167], 'lat', [-50, -40.8]);
    m_grid('box', 'fancy', 'tickdir', 'in', 'FontSize', 16);
    
    %
    hpc = m_pcolor(altdata.M2_Tasman.SUM.lon(linlonlims_aux), ...
                   altdata.M2_Tasman.SUM.lat(linlatlims_aux), ...
                   altdata.M2_Tasman.SUM.SSH(linlatlims_aux, linlonlims_aux) .*...
                                cos(altdata.M2_Tasman.SUM.PHS(linlatlims_aux, linlonlims_aux)));
    %
    hpc.LineStyle = 'none';


%
colormap(haxs_alt, clrmap.ssh);

%
caxis(haxs_alt, 15.*[-1, 1])

%
m_plot(mooringlocs.T1.longitude, mooringlocs.T1.latitude, '.k', 'MarkerSize', 60)

%
hcb_alt = colorbar(haxs_alt);
    hcb_alt.FontSize = 16;
    hcb_alt.Location = 'eastoutside';
    hcb_alt.Position = [0.845, 0.19, 0.0125, 0.5];
    hcb_alt.Label.String = '[mm]';
    hcb_alt.Label.FontSize = 26;
    hcb_alt.Label.Interpreter = 'Latex';
    hcb_alt.Label.Position(1) = 2.5;
    hcb_alt.Ticks = -15:5:15;

% Land
clr_land_aux = 0.5.*[1, 1, 1];
%
% % m_gshhs_l('patch', clr_land_aux);
m_gshhs_f('patch', clr_land_aux);

% % pause(1)

%
m_plot([lon_lims_plt(1), lon_lims_plt(2), lon_lims_plt(2), lon_lims_plt(1), lon_lims_plt(1)], ...
       [lat_lims_plt(2), lat_lims_plt(2), lat_lims_plt(1), lat_lims_plt(1), lat_lims_plt(2)], 'k')

% works when using  using 
hbla = get(haxs_alt, 'Children');
% % hbla(41).FaceColor = 0.8.*[1, 1, 1];
hbla(end).FaceColor = 0.8.*[1, 1, 1];

%
xlabel(haxs_alt, 'Longitude', 'Interpreter', 'Latex', 'FontSize', 18)
%
title(haxs_alt, 'M$_2$, mode-1 sea surface height', 'Interpreter', 'Latex', 'FontSize', 24)
        
        
% ---------------------------------------------
% ---------------------------------------------
%
axes(haxs_dsp)

    %
    m_proj('lambert', 'long', lon_lims_plt, 'lat', lat_lims_plt);
    m_grid('box', 'on', 'tickdir', 'in', 'FontSize', 16);
    hold on

    %
    hpc = m_pcolor(modelD.lon, modelD.lat, 1000*1025*modelD.Eps);
    hpc.LineStyle = 'none';

    %
    m_contour(modelD.lon, modelD.lat, modelD.Depth, [500:1000:5000], 'k')

%
colormap(haxs_dsp, clrmap.dissipation);

%
caxis(haxs_dsp, [0, 30])

%
hcb_dsp = colorbar(haxs_dsp);
    hcb_dsp.Position = [0.36, 0.2, 0.0125, 0.5];
    hcb_dsp.FontSize = 16;
    hcb_dsp.Label.String = '[mW m$^{-2}$]';
    hcb_dsp.Label.Interpreter = 'Latex';
    hcb_dsp.Label.FontSize = 18;

%
m_gshhs_f('patch', 0.5.*[1, 1, 1]);

%
lon_lims_small_aux = [148.9449, 149.06713];
lat_lims_small_aux = mooringlocs.T1.latitude + 0.035.*[-1, 1];

%
xlabel(haxs_dsp, 'Longitude', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel(haxs_dsp, 'Latitude', 'Interpreter', 'Latex', 'FontSize', 18)
%
title(haxs_dsp, ' \hspace{0.3cm} Dissipation', 'Interpreter', 'Latex', 'FontSize', 24)


%
hall_aux = get(haxs_dsp, 'Children');

%     %
%     hall_aux(16).Visible = 'off';
%     hall_aux(18).Visible = 'off';
%     hall_aux(19).Visible = 'off';
% %     %
hall_aux(61).Visible = 'off';
hall_aux(63).Visible = 'off';
hall_aux(64).Visible = 'off';

% -------------------------------------------------------

%
htxt_a = text(-0.012, 0.0262, 'a)');
    htxt_a.FontSize = 20;

%
htxt_b = text(0.035, 0.0262, 'b)');
    htxt_b.FontSize = 20;

%
% %     htxt_lbl = text(0.025, 0.027, {'\hspace{0.25cm}study';'\hspace{0.35cm}site';'(Fig. 2)'});
% %     htxt_lbl = text(0.03, 0.025, {'\hspace{0.25cm}study';'\hspace{0.375cm}site';'(Fig. 2)'});
htxt_lbl = text(0.025, 0.022, {'study';'site';'(Fig. 2)'});
    htxt_lbl.HorizontalAlignment = 'center';
    htxt_lbl.Interpreter = 'Latex';
    htxt_lbl.FontSize = 18;
        
        
% -------------------------------------------------------
%
haxs_annotation = axes('Position', [-0.5, -0.5, 2, 2]);
hold(haxs_annotation, 'on')
%
set(haxs_annotation, 'DataAspectRatio', [1, 1, 1])
set(haxs_annotation, 'XLim', [-1, 1], 'YLim', [-1, 1])

%
shaftwidth = 0.006;
arrwhdwidth = 0.012;
arrwhdlen = 0.03;
%
xcoordbase = [0, 0, (1-arrwhdlen), (1-arrwhdlen), 1, (1-arrwhdlen), (1-arrwhdlen), 0];
ycoordbase = [-shaftwidth/2, +shaftwidth/2, +shaftwidth/2, +shaftwidth/2 + arrwhdwidth, 0, (-shaftwidth/2 - arrwhdwidth), -shaftwidth/2, -shaftwidth/2];

%
x0 = [-0.2, -0.14, -0.001, 0.58];
y0 = [0.2, 0.2, 0.07, -0.1];
dir0 = [(-pi/2 - pi/3 - pi/60), (-pi/12)+(pi/60), (pi), (pi/2 + pi/3 - pi/60)];
mag0 = [0.2, 0.16, 0.272, 0.4];
% % %
% % x0 = [0, 0];
% % y0 = [0, -0.3];
% % dir0 = [pi/3, pi/3];
% % mag0 = [0.8, 0.5];
%
Narrows = length(x0);

%
scale_arrow = 1;

%
for i = 1:Narrows
    
    %
    if i==4
        %
        shaftwidth = 0.012;
        arrwhdwidth = 0.012;
        arrwhdlen = 0.03;
        %
        xcoordbase = [0, 0, (1-arrwhdlen), (1-arrwhdlen), 1, (1-arrwhdlen), (1-arrwhdlen), 0];
        ycoordbase = [-shaftwidth/2, +shaftwidth/2, +shaftwidth/2, +shaftwidth/2 + arrwhdwidth, 0, (-shaftwidth/2 - arrwhdwidth), -shaftwidth/2, -shaftwidth/2];

    end
    
    %
    x_arrow_aux = xcoordbase(:) .* scale_arrow * mag0(i);
	% Change neck position so it does not scale with magnitude
    x_arrow_aux([3:4, 6:7]) = x_arrow_aux(5) - arrwhdlen;
                                    
	%
    y_arrow_aux = ycoordbase(:);
% %     y_arrow_aux([]) = y_arrow_aux() .* scale_arrow * mag0(i);
    
    %
    rot_matrix = [cos(dir0(i)), -sin(dir0(i)); sin(dir0(i)), cos(dir0(i))];
    
    %
    xy_rot_aux = rot_matrix * [x_arrow_aux.'; y_arrow_aux.'];
    
    %
    x_arrow_plt = xy_rot_aux(1, :) + x0(i);
    y_arrow_plt = xy_rot_aux(2, :) + y0(i);
    
    %
    if i < 4
        clr_arrow = [0, 0, 0];
    else
        clr_arrow = [112, 191, 65]./255;
    end
    
    %
    hfill_aux = fill(haxs_annotation, x_arrow_plt, y_arrow_plt, clr_arrow);
    
end

%
htxtIT = text(0.325, 0.11, 'internal tide');
    htxtIT.Rotation = (180/pi)*dir0(4) - 180;
	htxtIT.FontSize = 24;

%
haxs_annotation.Visible = 'off';

        
%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure01.pdf'), 'Resolution', 300)

        
        
        
        
        
        
        
        
        
        
        
        
        
