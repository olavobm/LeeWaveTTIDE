%% Figure with NT1 and T1 data

clear
close all

%% Load data

% NT1
NT1data = load(fullfile(paper_directory(), 'data', 'shipboard', 'NT1', 'NT1_towyodata.mat'));
NT1data = NT1data.NT1data;

% T1 velocity
moorvel = load(fullfile(paper_directory(), 'data', 'mooring', 'level_2', 'T1_velocity_L2.mat'));
moorvel = moorvel.adcpdata;

% T1 dissipation
moorturb = load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_turbulence.mat'));
moorturb = moorturb.tchaindata;


%% Mooring locations

%
mooringlocs = load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));
mooringlocs = mooringlocs.mooringlocs;


%% Calculate cartesian (x, y) coordinates

%
lonlat0 = [149.0006, -41.3247];

%
[NT1data.bathy.x, ...
 NT1data.bathy.y] = lonlat2kmgrid(lonlat0, NT1data.bathy.lon, ...
                                           NT1data.bathy.lat);

%
for i = 1:NT1data.nrepeats
    
	%
    [NT1data.CTDdata.transects(i).x, ...
     NT1data.CTDdata.transects(i).y] = lonlat2kmgrid(lonlat0, ...
                                                     NT1data.CTDdata.transects(i).lon, ...
                                                     NT1data.CTDdata.transects(i).lat);
                                                 
    %
    [NT1data.LADCPdata.transects(i).x, ...
     NT1data.LADCPdata.transects(i).y] = lonlat2kmgrid(lonlat0, ...
                                                     NT1data.LADCPdata.transects(i).lon, ...
                                                     NT1data.LADCPdata.transects(i).lat);
                                                 
end

%
list_moor = {'T1', 'T2'};
%
for i = 1:length(list_moor)
	[mooringlocs.(list_moor{i}).x, ...
     mooringlocs.(list_moor{i}).y] = lonlat2kmgrid(lonlat0, ...
                                                   mooringlocs.(list_moor{i}).longitude, ...
                                                   mooringlocs.(list_moor{i}).latitude);
end




%% 
% -----------------------------
% --- GET AVERAGED VELOCITY ---
% ----- AND DO TIDAL FIT ------
% -----------------------------

%%

%
zlims_avg = [1450, 1800];

%
linzlimsavg = (NT1data.LADCPdata.z >= zlims_avg(1)) & ...
              (NT1data.LADCPdata.z <= zlims_avg(2));
          
%
timeavgvec = NaN(1, NT1data.nrepeats);
uavgvec = NaN(1, NT1data.nrepeats);

%
for i = 1:NT1data.nrepeats
    
    %
    timeavgvec(i) = nanmean(NT1data.LADCPdata.transects(i).time);
    
    %
    uavgvec(i) = nanmean(nanmean(NT1data.LADCPdata.transects(i).u(linzlimsavg, :), 1));
    
end

%
timeavgvec = timeavgvec - datenum(2015, 1, 1);


%%

%
mean_time = mean(timeavgvec);

%
x_var = timeavgvec(:) - mean_time;
%
y_var = uavgvec(:);

%
M2_freq = 24/12.42;

%
x_model_array = [ones(length(x_var), 1), x_var, ...
                 cos(2*pi*M2_freq.*x_var), sin(2*pi*M2_freq.*x_var)];

%
[mcoefs, mcoefs_bounds] = regress(y_var, x_model_array);


%% Compute a reconstruction to plot

%
time_grid = (18.45 : 0.01 : 19.5) - mean_time;

%
u_fit = mcoefs(1) + mcoefs(2).*time_grid(:) + ...
                    mcoefs(3).*cos(2*pi*M2_freq.*time_grid(:)) + ...
                    mcoefs(4).*sin(2*pi*M2_freq.*time_grid(:));


%
time_grid = time_grid + mean_time;


%% Compute potential density (relative to 2000 dbar,
% consistent with NT1) at T1

%
moorturb.data.sgth = gsw_sigma2(moorturb.data.SA, moorturb.data.CT);
%
moorturb.data.sgth = moorturb.data.sgth + 1000;


%% 
% -----------------------------
% - MEAN EPSILON IN TRANSECTS -
% -----------------------------

%
indplt = [1, 2];

%
epsi_1_aux = NT1data.CTDdata.transects(indplt(1)).procepsiOT(:);
epsi_2_aux = NT1data.CTDdata.transects(indplt(2)).procepsiOT(:);
%
epsibackground_aux = 1e-9;
%
epsi_1_aux(isnan(epsi_1_aux)) = epsibackground_aux;
epsi_2_aux(isnan(epsi_2_aux)) = epsibackground_aux;



    
%%
% -------------------------------------------------------
% --------------------- MAKE FIGURE ---------------------
% -------------------------------------------------------

%%

%
cmap.vel = load('colormap_velocity.mat');
cmap.epsilon = load('colormap_epsilon.mat');
%
cmap.vel = cmap.vel.cmapredblue;
cmap.epsilon = cmap.epsilon.cmaporange;

%
sgth_ctr_plt = 1036.7 : 0.0125 : 1036.95;


%%

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3032    0.2417    0.4227    0.5025];
%
haxs_NT1 = makeSubPlots(0.085, 0.4, 0.0125, ...
                        0.27, 0.12, 0.05, 2, 2);
haxs_timesr = axes('Position', [0.213, 0.85, 0.3, 0.1]);
haxs_T1 = makeSubPlots(0.765, 0.035, 0.0125, 0.27, 0.12, 0.05, 1, 2);
%
haxs_all = [haxs_NT1, haxs_timesr, haxs_T1];
hold(haxs_all, 'on')
   

    % ----------------------------------------------------
    %
    pcolor(haxs_NT1(1), NT1data.LADCPdata.transects(indplt(1)).x, ...
                        NT1data.LADCPdata.z, ...
                        NT1data.LADCPdata.transects(indplt(1)).u)
    %
    pcolor(haxs_NT1(2), NT1data.LADCPdata.transects(indplt(2)).x, ...
                        NT1data.LADCPdata.z, ...
                        NT1data.LADCPdata.transects(indplt(2)).u)
                   
	% ----------------------------------------------------
                    
    %
    zmatrix_1_aux = repmat(NT1data.CTDdata.z, 1, NT1data.CTDdata.transects(indplt(1)).nprofiles);
    zmatrix_2_aux = repmat(NT1data.CTDdata.z, 1, NT1data.CTDdata.transects(indplt(2)).nprofiles);
    
    % -------------------
    %
    epsi_aux = NT1data.CTDdata.transects(indplt(1)).procepsiOT(:);
    lok_aux = ~isnan(epsi_aux);
    %
    x_aux = NT1data.CTDdata.transects(indplt(1)).x(:);
    z_aux = zmatrix_1_aux(:);

    %
    scatter(haxs_NT1(3), x_aux(lok_aux), z_aux(lok_aux), 80, log10(epsi_aux(lok_aux)), 'filled')
   
    % -------------------
    %
    epsi_aux = NT1data.CTDdata.transects(indplt(2)).procepsiOT(:);
    lok_aux = ~isnan(epsi_aux);
    %
    x_aux = NT1data.CTDdata.transects(indplt(2)).x(:);
    z_aux = zmatrix_2_aux(:);

    %
    scatter(haxs_NT1(4), x_aux(lok_aux), z_aux(lok_aux), ...
                     80, log10(epsi_aux(lok_aux)), 'filled')

        
    % ----------------------------------------------------
        %
    % %         sgth_ref_bold = [1036.8, 1036.85];
    sgth_ref_bold = sgth_ctr_plt([11, 13]);

    %
    for i = [1, 3]
        %
        plot(haxs_NT1(i), NT1data.CTDdata.transects(indplt(1)).x(:), zmatrix_1_aux(:), ':k')
        %
        contour(haxs_NT1(i), NT1data.CTDdata.transects(indplt(1)).x, ...
                             zmatrix_1_aux, ...
                             NT1data.CTDdata.transects(indplt(1)).sgth, ...
                             sgth_ctr_plt, 'k', 'LineWidth', 1)
        %
        contour(haxs_NT1(i), NT1data.CTDdata.transects(indplt(1)).x, ...
                             zmatrix_1_aux, ...
                             NT1data.CTDdata.transects(indplt(1)).sgth, ...
                             sgth_ref_bold, 'k', 'LineWidth', 3)
    end

    %
    for i = [2, 4]
        %
        plot(haxs_NT1(i), NT1data.CTDdata.transects(indplt(2)).x(:), ...
                          zmatrix_2_aux(:), ':k')
        %
        contour(haxs_NT1(i), NT1data.CTDdata.transects(indplt(2)).x, ...
                             zmatrix_2_aux, ...
                             NT1data.CTDdata.transects(indplt(2)).sgth, ...
                             sgth_ctr_plt, 'k', 'LineWidth', 1)
        %
        contour(haxs_NT1(i), NT1data.CTDdata.transects(indplt(2)).x, ...
                             zmatrix_2_aux, ...
                             NT1data.CTDdata.transects(indplt(2)).sgth, ...
                             sgth_ref_bold, 'k', 'LineWidth', 3)
    end

    % ----------------------------------------------------
    % Plot T1 timeseries

    %
    time_plt_lims = [18.38, 18.94];
    %
    lvel_inlims = (moorvel.yday >= time_plt_lims(1)) & (moorvel.yday <= time_plt_lims(2));
    lsgth_inlims = (moorturb.data.yday >= time_plt_lims(1)) & (moorturb.data.yday <= time_plt_lims(2)); 
    lepsi_inlims = (moorturb.turbulence.yday >= time_plt_lims(1)) & (moorturb.turbulence.yday <= time_plt_lims(2)); 

    %
    pcolor(haxs_T1(1), moorvel.yday(lvel_inlims), moorvel.z, moorvel.u(:, lvel_inlims))
    pcolor(haxs_T1(2), moorturb.turbulence.yday(lepsi_inlims), moorturb.turbulence.zgrid, ...
                       log10(moorturb.turbulence.epsilon_gridded(:, lepsi_inlims)))
    %
    for i = 1:2
        %
        contour(haxs_T1(i), repmat(moorturb.data.yday(lsgth_inlims).', size(moorturb.data.depth, 1), 1), ...
                            moorturb.data.depth(:, lsgth_inlims), ...
                            moorturb.data.sgth(:, lsgth_inlims), sgth_ctr_plt, 'k')
        %
        contour(haxs_T1(i), repmat(moorturb.data.yday(lsgth_inlims).', size(moorturb.data.depth, 1), 1), ...
                            moorturb.data.depth(:, lsgth_inlims), ...
                            moorturb.data.sgth(:, lsgth_inlims), sgth_ref_bold, 'k', 'LineWidth', 2.5)
    end
                    
                    
%
for i = 1:length(haxs_all)
    shading(haxs_all(i), 'flat')
end

% ----------------------------------------------------
%
set(haxs_NT1, 'XLim', [-2.5, 1.8], 'YLim', [1400, 2025])
    %
    for i = 1:length(haxs_NT1)
        %
        plot(haxs_NT1(i), mooringlocs.T1.x.*[1, 1], [1400, 1950], '-k', 'LineWidth', 3);
        plot(haxs_NT1(i), mooringlocs.T2.x.*[1, 1], [1400, 1950], '-k', 'LineWidth', 3);
    end
    
    %
    for i = 1:length(haxs_NT1)
        axes(haxs_NT1(i))
            fill2Dbathy(NT1data.bathy.x, NT1data.bathy.depth);
    end
   
    %
    for i = 1:length(haxs_T1)
        axes(haxs_T1(i))
            fill2Dbathy(time_plt_lims, 1920.*[1, 1]);
    end

% ----------------------------------------------------
%
ship_time_1 = mean(NT1data.CTDdata.transects(indplt(1)).time(:, 3), 'omitnan') - datenum(2015, 1, 1);
ship_time_2 = mean(NT1data.CTDdata.transects(indplt(2)).time(:, end), 'omitnan') - datenum(2015, 1, 1);

    %
    for i = 1:2
        hp_mk_1 = plot(haxs_T1(i), ship_time_1, 1990, '^r', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
        hp_mk_2 = plot(haxs_T1(i), ship_time_2, 1990, '^g', 'MarkerSize', 16, 'MarkerFaceColor', 'g');
    end


% ----------------------------------------------------
% Add arrow for ship and some annotation
    
%
xtick = 0.25;
xtip = 0.15;
ytick = 15;
yneck = 15;
%
x0 = -2.15;
y0_1 = 2000;
x0_1 = x0;
%
xarrow_1 = [(x0_1 - xtick), (x0_1 + xtick), (x0_1 + xtick), (x0_1 + xtick + xtip), (x0_1 + xtick), (x0_1 + xtick), (x0_1 - xtick), (x0_1 - xtick)];
yarrow_1 = [(y0_1 - ytick/2), (y0_1 - ytick/2), (y0_1 - ytick/2 - yneck), ...
            y0_1, (y0_1 + ytick/2 + yneck), (y0_1 + ytick/2), (y0_1 + ytick/2), (y0_1 - ytick/2)];
%
for i = [1, 3]
    axes(haxs_NT1(i))
        fill(xarrow_1, yarrow_1, 'w')

    %
    htxt_aux = text(x0_1, y0_1-(3*yneck), 'ship');
        htxt_aux.FontSize = 14;
        htxt_aux.Color = [1, 1, 1];
        htxt_aux.HorizontalAlignment = 'center';

end

% ---------------
x0_2 = x0+xtip;
y0_2 = y0_1;
xarrow_2 = [(x0_2 + xtick), (x0_2 - xtick), (x0_2 - xtick), (x0_2 - xtick - xtip), (x0_2 - xtick), (x0_2 - xtick), (x0_2 + xtick), (x0_2 + xtick)];
yarrow_2 = [(y0_2 - ytick/2), (y0_2 - ytick/2), (y0_2 - ytick/2 - yneck), ...
            y0_2, (y0_2 + ytick/2 + yneck), (y0_2 + ytick/2), (y0_2 + ytick/2), (y0_2 - ytick/2)];
%
for i = [2, 4]
    axes(haxs_NT1(i))
        fill(xarrow_2, yarrow_2, 'w')

    %
    htxt_aux = text(x0_2, y0_2-(3*yneck), 'ship');
        htxt_aux.FontSize = 14;
        htxt_aux.Color = [1, 1, 1];
        htxt_aux.HorizontalAlignment = 'center';
end

%
axes(haxs_NT1(3))
    %
    htxt_aux = text(-0.95, 1980, ['$\langle \epsilon \rangle$ = ' num2str(mean(epsi_1_aux*1e8), '%.1f') '$\times 10^{-8}$']);
    htxt_aux.Interpreter = 'Latex';
    htxt_aux.Color = [1, 1, 1];
    htxt_aux.FontSize = 16;
%
axes(haxs_NT1(4))
    %
    htxt_aux = text(-0.95, 1980, ['$\langle \epsilon \rangle$ = ' num2str(mean(epsi_2_aux*1e8), '%.1f') '$\times 10^{-8}$']);
    htxt_aux.Interpreter = 'Latex';
    htxt_aux.Color = [1, 1, 1];
    htxt_aux.FontSize = 16;

% ----------------------------------------------------
% Set colormaps and color limits

%
colormap(haxs_T1(1), cmap.vel)
colormap(haxs_T1(2), cmap.epsilon)
%
set(haxs_T1(1), 'CLim', 0.3.*[-1, 1])
set(haxs_T1(2), 'CLim', [-9, -6])

%
colormap(haxs_NT1(1), cmap.vel)
colormap(haxs_NT1(2), cmap.vel)
colormap(haxs_NT1(3), cmap.epsilon)
colormap(haxs_NT1(4), cmap.epsilon)

%
set(haxs_NT1(1:2), 'CLim', 0.3.*[-1, 1])
set(haxs_NT1(3:4), 'CLim', [-9, -6])

% ----------------------------------------------------
 
%
xcb = 0.615;
wcb = 0.0125;

%
hcb_1 = colorbar(haxs_NT1(2));
    hcb_1.Position(1) = xcb;
    hcb_1.Position(3) = wcb;
    hcb_1.Ticks = -0.3:0.1:0.3;
    hcb_1.Label.Interpreter = 'Latex';
    hcb_1.Label.String = '[m s$^{-1}$]';
    hcb_1.FontSize = 12;
    hcb_1.Label.FontSize = 14;

%
hcb_2 = colorbar(haxs_NT1(4));
    hcb_2.Position(1) = xcb;
    hcb_2.Position(3) = wcb;
    hcb_2.Ticks = -9:1:-5;
    hcb_2.Label.Interpreter = 'Latex';
    hcb_2.Label.String = '$\log_{10}\epsilon$ [m$^2$ s$^{-3}$]';
    hcb_2.FontSize = 12;
    hcb_2.Label.FontSize = 14;

%
set(haxs_all, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
%
set([haxs_NT1, haxs_T1], 'YDir', 'reverse', 'YLim', [1400, 2025])
set(haxs_NT1, 'XLim', [-2.5, 1.8])
set(haxs_T1, 'XLim', time_plt_lims)

%
xaxsT1 = xlim(haxs_T1(1));
x_T1_ticks = 0 : 3 : 12;
x_T1_cell = cell(1, length(x_T1_ticks));
for i = 1:length(x_T1_ticks)
    x_T1_cell{i} = num2str(x_T1_ticks(i));
end
set(haxs_T1, 'XTick', xaxsT1(1) + x_T1_ticks./24, 'XTickLabel', x_T1_cell)
set(haxs_T1(1), 'XTickLabel', [])
%
set(haxs_NT1(1:2), 'XTickLabel', [])
set(haxs_NT1([2, 4]), 'YTickLabel', [])
set(haxs_all, 'Color', 0.8.*[1, 1, 1])

% ----------------------------------------------------

%
xlabel(haxs_NT1(3), '[km]', 'Interpreter', 'Latex', 'FontSize', 18)
xlabel(haxs_NT1(4), '[km]', 'Interpreter', 'Latex', 'FontSize', 18)
%
ylabel(haxs_NT1(1), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel(haxs_NT1(3), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 14)

%
xlabel(haxs_T1(2), {'Time [hours'; ['from yearday ' num2str(time_plt_lims(1)) ']']}, 'Interpreter', 'Latex', 'FontSize', 15)

% ----------------------------------------------------
% Now add a plot with the timeseries and fit
%
axes(haxs_timesr)

    %
    NT1_1_min = min(NT1data.CTDdata.transects(indplt(1)).time(:)) - datenum(2015, 1, 1);
    NT1_1_max = max(NT1data.CTDdata.transects(indplt(1)).time(:)) - datenum(2015, 1, 1);
    %
    NT1_2_min = min(NT1data.CTDdata.transects(indplt(2)).time(:)) - datenum(2015, 1, 1);
    NT1_2_max = max(NT1data.CTDdata.transects(indplt(2)).time(:)) - datenum(2015, 1, 1);

    %
    xlim(haxs_timesr, [18.44, 19.54]);
    ylim(haxs_timesr, [-0.125, 0.1]);
    sliceShading([NT1_1_min, NT1_1_max; ...
                  NT1_2_min, NT1_2_max], ...
                  [1, 0, 0; 0, 1, 0], 1);

    %
    plot(haxs_timesr, timeavgvec, uavgvec, '.-b', 'MarkerSize', 30)
    plot(haxs_timesr, time_grid, u_fit, '--b')

    %
    plot(haxs_timesr, [18.44, 19.54], [0, 0], '-k')
    plot(haxs_timesr, [18.44, 19.54], mcoefs(1).*[1, 1], '--b')

%
xlim(haxs_timesr, [18.44, 19.54]);
ylim(haxs_timesr, [-0.105, 0.05]);
set(haxs_timesr, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
h_xlbl = xlabel(haxs_timesr, 'Time [yearday]', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel(haxs_timesr, '[m s$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', 14)
%
h_xlbl.Position(2) = -0.16;

% ----------------------------------------------------
% Add arrows


% ------------------------------------------------------------
% Plot arrows for the mean flow

%
shaftwidth = 0.04;
arrwhdwidth = 0.06;
arrwhdlen = 0.16;
%
xcoordbase = [0, 0, (1-arrwhdlen), (1-arrwhdlen), 1, (1-arrwhdlen), (1-arrwhdlen), 0];
ycoordbase = [-shaftwidth/2, +shaftwidth/2, +shaftwidth/2, +shaftwidth/2 + arrwhdwidth, 0, (-shaftwidth/2 - arrwhdwidth), -shaftwidth/2, -shaftwidth/2];

% --------------------------------
%
axslength = 0.35;
axs_arrw_1 = axes('Position', [0.05, 0.7, axslength, axslength]);
%
xy_coordlims = [xcoordbase(:).'; ycoordbase(:).'];
%
rot_ang = 240;
rot_matrix = [cosd(rot_ang), -sind(rot_ang); sind(rot_ang), cosd(rot_ang)];
%
xy_rotcoord = rot_matrix * xy_coordlims;

%
hfill_1 = fill(axs_arrw_1, xy_rotcoord(1, :), xy_rotcoord(2, :), 'r');
hfill_1.LineStyle = 'none';

%
set(axs_arrw_1, 'DataAspectRatio', [1, 1, 1])
set(axs_arrw_1, 'XLim', [-1, 1], 'YLim', [-1, 1])

%
axs_arrw_1.Visible = 'off';

% --------------------------------
%
axslength = 0.4;
axs_arrw_2 = axes('Position', [0.062, 0.65, axslength, axslength]);
%
xy_coordlims = [xcoordbase(:).'; ycoordbase(:).'];
%
rot_ang = -45;
rot_matrix = [cosd(rot_ang), -sind(rot_ang); sind(rot_ang), cosd(rot_ang)];
%
xy_rotcoord = rot_matrix * xy_coordlims;

%
hfill_2 = fill(axs_arrw_2, xy_rotcoord(1, :), xy_rotcoord(2, :), 'g');
hfill_2.LineStyle = 'none';

%
set(axs_arrw_2, 'DataAspectRatio', [1, 1, 1])
set(axs_arrw_2, 'XLim', [-1, 1], 'YLim', [-1, 1])

%
axs_arrw_2.Visible = 'off';



% ----------------------------------------------------
% Add letter labelling to subplots

%
text(haxs_timesr, 18.45, 0.085, 'a)', 'Interpreter', 'Latex', 'FontSize', 16, 'HorizontalAlignment', 'left')
text(haxs_timesr, 18.55, 0.085, 'Cross-shore velocity ($U$)', 'Interpreter', 'Latex', 'FontSize', 16, 'HorizontalAlignment', 'left')

%
text(haxs_NT1(1), -2.4, 1350, 'b) $U$', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_NT1(2), -2.4, 1350, 'c) $U$', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_NT1(3), -2.4, 1350, 'e)', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_NT1(4), -2.4, 1350, 'f)', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_NT1(3), -1.9, 1350, '$\epsilon$', 'Interpreter', 'Latex', 'FontSize', 26, 'HorizontalAlignment', 'left')
text(haxs_NT1(4), -1.9, 1350, '$\epsilon$', 'Interpreter', 'Latex', 'FontSize', 26, 'HorizontalAlignment', 'left')

%
text(haxs_T1(1), 18.385, 1350, 'd) $U$ at T1', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_T1(2), 18.385, 1350, 'g) $\hphantom{\epsilon}$ at T1', 'Interpreter', 'Latex', 'FontSize', 19, 'HorizontalAlignment', 'left')
text(haxs_T1(2), 18.465, 1350, '$\epsilon$', 'Interpreter', 'Latex', 'FontSize', 26, 'HorizontalAlignment', 'left')



%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure06.pdf'), 'Resolution', 300)

