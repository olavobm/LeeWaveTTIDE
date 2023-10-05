%%
clear
close all

%%
% -----------------------------
% --------- LOAD DATA ---------
% -----------------------------

%
T1ellison = load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_turbulence.mat'));
T1ellison = T1ellison.tchaindata;

%%

%
T1thorpe = load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_turbulence_Thorpe.mat'));
T1thorpe = T1thorpe.tchaindata;


%%
% -----------------------------
% ------- FILL NaNs IN --------
% -- OVERTURN-BASED EPSILON ---
% -----------------------------

%
epsi_background_fill = 1e-9;

%
T1thorpe.epsilon_original = T1thorpe.epsilon;
%
T1thorpe.epsilon(isnan(T1thorpe.epsilon)) = epsi_background_fill;


%%
% -----------------------------
% --- SIMPLE DEPTH AVERAGING --
% -----------------------------

%%

%
datagridded.zlimsavg = [1710, 1900];


%% Ellison

%
datagridded.Ellison.dt_data = T1ellison.turbulence.dt;
datagridded.Ellison.dtime_data = T1ellison.turbulence.dtime;
datagridded.Ellison.yday_data = T1ellison.turbulence.yday;

%
lin_zlims = (T1ellison.turbulence.zgrid >= datagridded.zlimsavg(1)) & ...
            (T1ellison.turbulence.zgrid <= datagridded.zlimsavg(2));

%
datagridded.Ellison.epsilon_depthavg = mean(T1ellison.turbulence.epsilon_gridded(lin_zlims, :), 1);


%% Thorpe

%
datagridded.Thorpe.dt_data = T1thorpe.dt;
datagridded.Thorpe.dtime_data = T1thorpe.dtime;
datagridded.Thorpe.yday_data = T1thorpe.yday;

%
lin_zlims = (T1thorpe.depth >= datagridded.zlimsavg(1)) & ...
            (T1thorpe.depth <= datagridded.zlimsavg(2));

%
datagridded.Thorpe.epsilon_depthavg = mean(T1thorpe.epsilon(lin_zlims, :), 1);


%%
% -----------------------------
% ----- SIMPLE SEMIDIURNAL ----
% --------- AVERAGING ---------
% -----------------------------


%%

%
datagridded.dt = days(0.5);
%
datagridded.dtime = datetime(2015, 01, 19) : datagridded.dt : datetime(2015, 03, 03);
datagridded.dtime.TimeZone = 'America/Los_Angeles';
%
datagridded.yday = datenum(datagridded.dtime) - datenum(2015, 1, 1);
% % datagridded.dtime = 18 : datagridded.dtlow : 61;


%% Ellison

%
ind_start = find(datagridded.Ellison.dtime_data >= (datagridded.dtime(1) - (datagridded.dt/2)), 1, 'first');
ind_stop = find(datagridded.Ellison.dtime_data < (datagridded.dtime(end) + (datagridded.dt/2)), 1, 'last');
%
ind_get = ind_start : 1 : ind_stop;

%
Nrows = datagridded.dt/datagridded.Ellison.dt_data;
Ncols = length(ind_get)/Nrows;
    
%
epsi_aux = reshape(datagridded.Ellison.epsilon_depthavg(ind_get), Nrows, Ncols);

%
datagridded.Ellison.epsilon = mean(epsi_aux, 1);
   

%% Thorpe

%
ind_start = find(datagridded.Thorpe.dtime_data >= (datagridded.dtime(1) - (datagridded.dt/2)), 1, 'first');
ind_stop = find(datagridded.Thorpe.dtime_data < (datagridded.dtime(end) + (datagridded.dt/2)), 1, 'last');
%
ind_get = ind_start : 1 : ind_stop;

%
Nrows = datagridded.dt/datagridded.Thorpe.dt_data;
Ncols = length(ind_get)/Nrows;
    
%
epsi_aux = reshape(datagridded.Thorpe.epsilon_depthavg(ind_get), Nrows, Ncols);

%
datagridded.Thorpe.epsilon = mean(epsi_aux, 1);


%%
% -----------------------------
% ---- COMPUTE REGRESSION -----
% ------ BETWEEN METHODS ------
% -----------------------------

%
xvar = log10(datagridded.Thorpe.epsilon(:));
yvar = log10(datagridded.Ellison.epsilon(:));

%
[mcoef, mcoef_int] = regress(yvar, [ones(length(xvar), 1), xvar]);




%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ------------------------- MAKE FIGURE -------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


%%

%
yday_lims_data = [38.85, 40.1];

%
linlims_Thorpe_data = (T1thorpe.yday >= yday_lims_data(1)) & ...
                      (T1thorpe.yday <= yday_lims_data(2));

%
linlims_Ellison_data = (T1ellison.data.yday >= yday_lims_data(1)) & ...
                       (T1ellison.data.yday <= yday_lims_data(2));
%
linlims_Ellison_turb = (T1ellison.turbulence.yday >= yday_lims_data(1)) & ...
                       (T1ellison.turbulence.yday <= yday_lims_data(2));


%%

T1ellison.data.sigma2 = gsw_sigma2(T1ellison.data.SA, T1ellison.data.CT);


%%

%
axs_scatter_lims = [8e-10, 5e-8];

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.2909    0.4125    0.4341    0.3125];
%
haxs_pcolors = makeSubPlots(0.08, 0.42, 0.055, 0.15, 0.2, 0.02, 2, 2);
haxs_scatter = axes('Position', [0.535, 0.2, 0.65, 0.65]);
%
haxs_all = [haxs_pcolors, haxs_scatter];
%
hold(haxs_pcolors, 'on')
hold(haxs_scatter, 'on')

    % ------------------------------------------------
    %
    pcolor(haxs_pcolors(1), T1thorpe.yday(linlims_Thorpe_data), T1thorpe.depth, T1thorpe.LT(:, linlims_Thorpe_data))
    pcolor(haxs_pcolors(3), T1thorpe.yday(linlims_Thorpe_data), T1thorpe.depth, log10(T1thorpe.epsilon_original(:, linlims_Thorpe_data)))
    
    %
    pcolor(haxs_pcolors(2), T1ellison.turbulence.yday(:, linlims_Ellison_turb), T1ellison.turbulence.depth(:, linlims_Ellison_turb), T1ellison.turbulence.LE(:, linlims_Ellison_turb))
    pcolor(haxs_pcolors(4), T1ellison.turbulence.yday(:, linlims_Ellison_turb), T1ellison.turbulence.zgrid, log10(T1ellison.turbulence.epsilon_gridded(:, linlims_Ellison_turb)))

    %
    for i = 1:length(haxs_pcolors)
        shading(haxs_pcolors(i), 'flat')
    end
    
    % Depth of last thermistor where spacing becomes bigger
    z_edgetop = 1710;
    %
    plot(haxs_pcolors(1), yday_lims_data, z_edgetop.*[1, 1], '-r', 'LineWidth', 2)
    plot(haxs_pcolors(3), yday_lims_data, z_edgetop.*[1, 1], '-r', 'LineWidth', 2)
    
    %
    sigma2_ctrs = 36.6 : 0.025 : 36.95;
    %
    for i = 1:length(haxs_pcolors)
        contour(haxs_pcolors(i), repmat(T1ellison.data.yday(linlims_Ellison_data).', length(T1ellison.SN), 1), ...
                                 T1ellison.data.depth(:, linlims_Ellison_data), ...
                                 T1ellison.data.sigma2(:, linlims_Ellison_data), ...
                                 sigma2_ctrs, '-k', 'LineWidth', 1)
    end
    
    % ------------------------------------------------
    %
    plot(haxs_scatter, datagridded.Thorpe.epsilon, ...
                       datagridded.Ellison.epsilon, '.k', 'MarkerSize', 24)
	%
    plot(haxs_scatter, axs_scatter_lims, axs_scatter_lims, '--k', 'LineWidth', 2)
    plot(haxs_scatter, axs_scatter_lims, (10^mcoef(1))*axs_scatter_lims.^mcoef(2), '-r', 'LineWidth', 2)
    
% ------------------------------------------------
%
xdb_1 = 0.3075;
xdb_2 = 0.585;
wdb = 0.015;
%
hcb_1 = colorbar(haxs_pcolors(1));
hcb_2 = colorbar(haxs_pcolors(3));
%
hcb_3 = colorbar(haxs_pcolors(2));
hcb_4 = colorbar(haxs_pcolors(4));
%
hcb_1.Position(1) = xdb_1;
hcb_2.Position(1) = xdb_1;
%
hcb_3.Position(1) = xdb_2;
hcb_4.Position(1) = xdb_2;
%
hcb_3.Position(3) = wdb;
hcb_4.Position(3) = wdb;
%
%
hcb_3.Label.Interpreter = 'Latex';
hcb_4.Label.Interpreter = 'Latex';
%
hcb_3.Label.String = '[m]';
hcb_4.Label.String = '$\log_{10} \epsilon~[\mathrm{m}^2~\mathrm{s}^{-3}]$';

%
cmapepsi = load('colormap_epsilon.mat');
cmapepsi = cmapepsi.cmaporange;
colormap(haxs_pcolors(3), cmapepsi)
colormap(haxs_pcolors(4), cmapepsi)

% ------------------------------------------------
    
%
set(haxs_all, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'YTick', 1500:200:1900)
%
set(haxs_pcolors, 'FontSize', 12, 'YDir', 'reverse', 'XLim', [38.9345, 40.0345], 'YLim', [1400, 1920])
set(haxs_pcolors(1:2), 'XTickLabel', [])
set(haxs_pcolors([2, 4]), 'YTickLabel', [])
set(haxs_pcolors(1), 'CLim', [0, 50])
set(haxs_pcolors(3), 'CLim', [-9, -6]);
set(haxs_pcolors(2), 'CLim', [0, 25])
set(haxs_pcolors(4), 'CLim', [-9, -7]);
%
set(haxs_pcolors, 'Color', 0.5.*[1, 1, 1])
%
set(haxs_scatter, 'DataAspectRatio', [1, 1, 1])
set(haxs_scatter, 'XScale', 'log', 'YScale', 'log', 'XLim', axs_scatter_lims, 'YLim', axs_scatter_lims)
set(haxs_scatter, 'XTick', 10.^(-10:1:-6), 'YTick', 10.^(-10:1:-6))
    
% ------------------------------------------------
%
ylabel(haxs_pcolors(1), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel(haxs_pcolors(3), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 14)
%
xlabel(haxs_pcolors(3), 'Time [yearday]', 'Interpreter', 'Latex', 'FontSize', 14)
xlabel(haxs_pcolors(4), 'Time [yearday]', 'Interpreter', 'Latex', 'FontSize', 14)

%
title(haxs_pcolors(1), 'Thorpe-based', 'Interpreter', 'Latex', 'FontSize', 18)
title(haxs_pcolors(2), 'Ellison-based', 'Interpreter', 'Latex', 'FontSize', 18)
%
xlabel(haxs_scatter, '$\epsilon$ from $L_T$ [m$^2$ s$^{-3}$]', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel(haxs_scatter, '$\epsilon$ from $L_E$ [m$^2$ s$^{-3}$]', 'Interpreter', 'Latex', 'FontSize', 16)
title(haxs_scatter, 'Semidiurnal, depth-averages', 'Interpreter', 'Latex', 'FontSize', 14)

% ------------------------------------------------
% Add letter labelling

% ------------------
%
xrectcoord = [38.94, 39.225, 39.225, 38.94, 38.94];
yrectcoord = [1406, 1406, 1530, 1530, 1406];
%
for i = 1:length(haxs_pcolors)
    fill(haxs_pcolors(i), xrectcoord, yrectcoord, 'w')
end

%
textlblFS = 16;

%
letters_aux = {'a)', 'b)', 'c)', 'd)'};
%
for i = 1:length(haxs_pcolors)
	%
    text(haxs_pcolors(i), 38.96, mean(yrectcoord(2:3)), letters_aux{i}, ...
                          'FontSize', textlblFS, 'HorizontalAlignment', 'left')
end
text(haxs_pcolors(1), 39.1-0.032, mean(yrectcoord(2:3))+8, '$L_T$', 'FontSize', textlblFS+2, 'Interpreter', 'Latex', 'HorizontalAlignment', 'left')
text(haxs_pcolors(2), 39.1-0.032, mean(yrectcoord(2:3))+8, '$L_E$', 'FontSize', textlblFS+2, 'Interpreter', 'Latex', 'HorizontalAlignment', 'left')
text(haxs_pcolors(3), 39.1, mean(yrectcoord(2:3)), '$\epsilon$', 'FontSize', textlblFS+12, 'Interpreter', 'Latex', 'HorizontalAlignment', 'left')
text(haxs_pcolors(4), 39.1, mean(yrectcoord(2:3)), '$\epsilon$', 'FontSize', textlblFS+12, 'Interpreter', 'Latex', 'HorizontalAlignment', 'left')

% ------------------
%
xrectcoord = [8.25e-10, 8.25e-10, 1.4e-9, 1.4e-9, 8.25e-10];
yrectcoord = [3e-8, 4.8e-8, 4.8e-8, 3e-8, 3e-8];
%
fill(haxs_scatter, xrectcoord, yrectcoord, 'w')
%
text(haxs_scatter, 9e-10, 3.8e-8, 'e)', 'FontSize', textlblFS, 'HorizontalAlignment', 'left')


%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure03.pdf'), 'Resolution', 300)




