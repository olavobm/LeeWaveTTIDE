%% Figure with Fast-CTD transect and T1 data

clear
close all


%% Load data

% FastCTD
FCTDdata = load(fullfile(paper_directory(), 'data', 'shipboard', 'FastCTD', 'FastCTD_towyo.mat'));
FCTDdata = FCTDdata.FCTDtransect;


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


%% Compute potential density relative to 2000 dbar

% Mooring T1
%
moorturb.data.sgth = gsw_sigma2(moorturb.data.SA, moorturb.data.CT);
%
moorturb.data.sgth = moorturb.data.sgth + 1000;


%% Compute distance of T1 correspondent to the FCTD transect

%
lon_ref = 149.00300139959;

% In km
mooringlocs.T1.xdist = distAlongTransect([FCTDdata.bathy.lon(1), FCTDdata.zgrid.lat0; ...
                                          FCTDdata.bathy.lon(end), FCTDdata.zgrid.lat0], ...
                                         [mooringlocs.T1.longitude, mooringlocs.T1.latitude], true);
x_bump = distAlongTransect([FCTDdata.bathy.lon(1), FCTDdata.zgrid.lat0; ...
                            FCTDdata.bathy.lon(end), FCTDdata.zgrid.lat0], ...
                           [lon_ref, FCTDdata.zgrid.lat0], true);
%
mooringlocs.T1.xdist = mooringlocs.T1.xdist - x_bump;
                       

%% 
% -----------------------------
% --- COMPUTE INTERNAL-WAVE ---
% ------ CHARACTERISTIC -------
% -----------------------------

%% Compute N2 profile

%
timeavg_lims = [41, 48];

%
linavglims = (moorturb.data.yday >= timeavg_lims(1)) & ...
             (moorturb.data.yday <= timeavg_lims(2));
         
%
SA_mean = mean(moorturb.data.SA(:, linavglims), 2, 'omitnan');
CT_mean = mean(moorturb.data.CT(:, linavglims), 2, 'omitnan');
p_mean = mean(moorturb.data.pressure(:, linavglims), 2, 'omitnan');

%
[N2_mean, ~] = gsw_Nsquared(SA_mean, CT_mean, p_mean, mooringlocs.T1.latitude);

%
z_mean = mean(moorturb.data.depthmid(:, linavglims), 2, 'omitnan');


%% 

%
x_array = [ones(length(z_mean), 1), z_mean(:)];

%
mcoefs = regress(N2_mean, x_array);

%
iwchar.freq = (2*pi)*(24/12.42)./(24*3600);

%
iwchar.zgrid = 1000:1:2100;

%
iwchar.N2prof = mcoefs(1) + mcoefs(2).*iwchar.zgrid;

%
iwchar.lat0 = mooringlocs.T1.latitude;
iwchar.f0 = gsw_f(iwchar.lat0)^2;

%
iwchar.charslope = sqrt((iwchar.freq^2 - iwchar.f0^2) ./ ...
                        (iwchar.N2prof - iwchar.freq^2));
      
                    
%% Trace ray

%
iwchar.ray.xz0 = [0, 1800];

%
iwchar.ray.Nsteps = 400;
iwchar.ray.dxstep = 10;

%
iwchar.ray.x = NaN(iwchar.ray.Nsteps, 1);
iwchar.ray.z = NaN(iwchar.ray.Nsteps, 1);
%
iwchar.ray.x(1) = iwchar.ray.xz0(1);
iwchar.ray.z(1) = iwchar.ray.xz0(2);
    
%
for i = 2:iwchar.ray.Nsteps

    %
    angTan_aux = interp1(iwchar.zgrid, iwchar.charslope, iwchar.ray.z(i-1));

    %
    iwchar.ray.x(i) = iwchar.ray.x(i-1) - iwchar.ray.dxstep;
    iwchar.ray.z(i) = iwchar.ray.z(i-1) - (iwchar.ray.dxstep * angTan_aux);
end

%
iwchar.ray.x = iwchar.ray.x ./ 1000;


          
    
%%
% -------------------------------------------------------
% --------------------- MAKE FIGURE ---------------------
% -------------------------------------------------------


%%  
% ------------------------------------------------------------
%
% time_lims = [42.65, 42.83];
time_lims = [42.6, 43];

%
xtw_lims = [-3.5, 2];

% ------------------------------------------------------------
%
linT1vel_aux = (moorvel.yday >= time_lims(1)) & ...
               (moorvel.yday <= time_lims(2));
%
linT1epsi_aux = (moorturb.turbulence.yday >= time_lims(1)) & ...
                (moorturb.turbulence.yday <= time_lims(2));
%
linT1sgth_aux = (moorturb.data.yday >= time_lims(1)) & ...
                (moorturb.data.yday <= time_lims(2));


% Towyo transect to plot
indtransect_plt = 2;

% Time limits for the towyo transect
ind_pos_ref = [1, 38, ...
               (size(FCTDdata.zgrid.transects(indtransect_plt).time, 2) - 5)];
twtimerefs = nanmean(FCTDdata.zgrid.transects(indtransect_plt).time(:, ind_pos_ref), 1);
xpostwrefs = nanmean(FCTDdata.zgrid.transects(indtransect_plt).x(:, ind_pos_ref), 1);


%%

%
cmap.vel = load('colormap_velocity.mat');
cmap.epsilon = load('colormap_epsilon.mat');
%
cmap.vel = cmap.vel.cmapredblue;
cmap.epsilon = cmap.epsilon.cmaporange;


%%

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3200    0.1083    0.3695    0.5308];
%
haxs_moor = makeSubPlots(0.6, 0.15, 0.1, ...
                         0.1, 0.15, 0.02, 1, 3);
haxs_FCTD = makeSubPlots(0.125, 0.6125, 0.1, ...
                         0.35, 0.35, 0.02, 1, 1);
%
haxs_all = [haxs_moor, haxs_FCTD];
hold(haxs_all, 'on')
  

    % -----------------------------------------------------  
    %
    indtransect_plt = 2;

    %
    zmatrix_aux = repmat(FCTDdata.zgrid.z, 1, FCTDdata.zgrid.profsperrepeat(indtransect_plt));
    %
    x_FCTD_vec = FCTDdata.zgrid.transects(indtransect_plt).x(:);
    y_FCTD_vec = zmatrix_aux(:);
    z_FCTD_vec = log10(FCTDdata.zgrid.transects(indtransect_plt).OTs.epsiOTflagged(:));
    %
    ldataOK_aux = ~isnan(z_FCTD_vec);
    %

    scatter(haxs_FCTD, x_FCTD_vec(ldataOK_aux), ...
                       y_FCTD_vec(ldataOK_aux), 30,...
                       z_FCTD_vec(ldataOK_aux), 'filled')

    %
    contour(haxs_FCTD, FCTDdata.zgrid.transects(indtransect_plt).x, ...
                       repmat(FCTDdata.zgrid.z, 1, FCTDdata.zgrid.profsperrepeat(indtransect_plt)), ...
                       FCTDdata.zgrid.transects(indtransect_plt).sigma2, 20, 'k')

    %
    plot(haxs_FCTD, mooringlocs.T1.xdist.*[1, 1], [1200, 2000], '-k', 'LineWidth', 2)

    %
    pcolor(haxs_moor(2), moorvel.yday(linT1vel_aux), moorvel.z, 100*moorvel.u(:, linT1vel_aux))
    pcolor(haxs_moor(3), moorvel.yday(linT1vel_aux), moorvel.z, 100*moorvel.v(:, linT1vel_aux))

    %
    pcolor(haxs_moor(1), moorturb.turbulence.yday(linT1epsi_aux), ...
                         moorturb.turbulence.zgrid, ...
                         log10(moorturb.turbulence.epsilon_gridded(:, linT1epsi_aux)))

    %
    sgthctr_lvls = 1036.6 : 0.0125 : 1036.95;
    %
    for i = 1:length(haxs_moor)
            %
            contour(haxs_moor(i), repmat(moorturb.data.yday(linT1sgth_aux).', size(moorturb.data.depth, 1), 1), ...
                                  moorturb.data.depth(:, linT1sgth_aux), ...
                                  moorturb.data.sgth(:, linT1sgth_aux), sgthctr_lvls, 'k')
    end
     
    %
    for i2 = 1
        plot(haxs_FCTD, iwchar.ray.x, iwchar.ray.z, '--k', 'LineWidth', 5)
    end

%
for i = 1:length(haxs_all)
    shading(haxs_all(i), 'flat')
end

                           
% -----------------------------------------------------

%
colormap(haxs_FCTD, cmap.epsilon)
colormap(haxs_moor(1), cmap.epsilon)
%
colormap(haxs_moor(2), cmap.vel)
colormap(haxs_moor(3), cmap.vel)

%
axes(haxs_FCTD)
    fill2Dbathy(FCTDdata.bathy.x, FCTDdata.bathy.depth);

%
xcb_1 = 0.4;
xcb_2 = 0.86;
wcb = 0.016;
%
hcb_1 = colorbar(haxs_FCTD);
    hcb_1.FontSize = 14;
    hcb_1.Position(1) = xcb_1;
    hcb_1.Position(3) = wcb;
    hcb_1.Label.Interpreter = 'Latex';
    hcb_1.Label.String = '$\log_{10}\epsilon~[\mathrm{m}^2~\mathrm{s}^{-3}]$';
    hcb_1.Label.FontSize = 18;
    hcb_1.Ticks = -9:1:-6;
%
hcb_2 = colorbar(haxs_moor(1));
    hcb_2.FontSize = 14;
    hcb_2.Position(1) = xcb_2;
    hcb_2.Position(3) = wcb;
    hcb_2.Label.Interpreter = 'Latex';
    hcb_2.Label.String = '$\log_{10}\epsilon~[\mathrm{m}^2~\mathrm{s}^{-3}]$';
    hcb_2.Label.FontSize = 16;
    hcb_2.Ticks = -9:1:-6;
%
hcb_3 = colorbar(haxs_moor(3));
    hcb_3.FontSize = 14;
    hcb_3.Position(1) = xcb_2;
    hcb_3.Position(3) = wcb;
    hcb_3.Position(4) = sum(haxs_moor(2).Position([2, 4])) - haxs_moor(3).Position(2);
    hcb_3.Label.Interpreter = 'Latex';
    hcb_3.Label.String = '[m s$^{-1}$]';
    hcb_3.Label.FontSize = 20;
    hcb_3.Ticks = -20:10:20;
    hcb_3.TickLabels = {'-0.2', '-0.1', '0', '0.1', '0.2'};
       
%
set([haxs_FCTD, haxs_moor(1)], 'CLim', [-9, -6])
set(haxs_moor(2:3), 'CLim', 22.*[-1, 1])
    

% -----------------------------------------------------

%
set(haxs_moor(1:2), 'XTickLabel', [])
    
%
set(haxs_all, 'FontSize', 14, 'Box', 'on', 'YDir', 'reverse', 'YLim', [1400, 1920])
set(haxs_moor, 'XLim', time_lims)
set([haxs_FCTD, haxs_moor], 'YTick', 1500:200:1900)
set(haxs_FCTD, 'XTick', -5:1:3)
set(haxs_FCTD, 'XLim', xtw_lims)
set(haxs_all, 'Color', 0.8.*[1, 1, 1])
%
set(haxs_moor, 'XTick', time_lims(1) + [0:3:9]./24)
set(haxs_moor(3), 'XTickLabel', {'0', '3', '6', '9'})
    

% -----------------------------------------------------
%
xlabel(haxs_moor(3), {'Time [hours'; ['from yearday ' num2str(haxs_moor(1).XLim(1)) ']']}, 'Interpreter', 'Latex', 'FontSize', 17)
%    
xlabel(haxs_FCTD, '[km]', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel(haxs_FCTD, 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', 14)
       
    
% ------------------------------------------------------------
% Add red arrows to reference ship location

%
haxs_axs_FCTD = axes('Position', [haxs_FCTD.Position(1), 0.64, haxs_FCTD.Position(3), 0.075]);
haxs_axs_moor = axes('Position', [haxs_moor(1).Position(1), 0.89, haxs_moor(1).Position(3), 0.075]);
hold([haxs_axs_FCTD, haxs_axs_moor], 'on')

%
nmarkref = length(xpostwrefs);
markerSZ = 16;
y0markersSZref = 1760;
%
plot(haxs_axs_FCTD, xpostwrefs, y0markersSZref.*ones(1, nmarkref), 'vk', ...
                      'MarkerSize', markerSZ, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
plot(haxs_axs_moor, twtimerefs, y0markersSZref.*ones(1, nmarkref), 'vk', ...
                      'MarkerSize', markerSZ, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')

%
axes(haxs_axs_FCTD)
    %
    x0topback = -0.5;
    xlen2neck = 1.5;
    xlenhead = 0.4;
    xarrowvec = [x0topback, ...
                 (x0topback + xlen2neck), (x0topback + xlen2neck), ...
                 (x0topback + xlen2neck + xlenhead), ...
                 (x0topback + xlen2neck), (x0topback + xlen2neck), ...
                 x0topback, x0topback];

    %
    y0topback = 1700;
    y0topback = y0markersSZref;
    y0topback = y0markersSZref - 120;
    %
    ywdtip = 100;
    ywdneck = 120;
    yarrowvec = [y0topback, y0topback, ...
                 (y0topback - ywdneck), ...
                 (y0topback + ywdtip), ...
                 (y0topback + 2*ywdtip + ywdneck), ...
                 (y0topback + 2*ywdtip), (y0topback + 2*ywdtip), y0topback];
    %
    hfill_aux = fill(xarrowvec, yarrowvec,'r');
        hfill_aux.LineStyle = 'none';

    %
    htxt_SHIP = text(-0.4, 1735, 'S H I P');
        htxt_SHIP.FontSize = 14;

%
xlim(haxs_axs_FCTD, xtw_lims)
xlim(haxs_axs_moor, time_lims)
%
set([haxs_axs_FCTD, haxs_axs_moor], 'Visible', 'off')
    

% ------------------------------------------------------
% Add titles
%
titleFS = 22;
% titley0ref = 1290;
titley0ref = 1270;

%
htxt1_aux = text(haxs_FCTD, mean(xtw_lims), titley0ref, 'Fast-CTD towyo');
%
htxt2_aux = text(haxs_moor(1), mean(time_lims), titley0ref, 'T1 mooring');

%
htxt1_aux.HorizontalAlignment = 'center';
htxt2_aux.HorizontalAlignment = 'center';
%
htxt1_aux.Interpreter = 'Latex';
htxt2_aux.Interpreter = 'Latex';
%
htxt1_aux.FontSize = titleFS;
htxt2_aux.FontSize = titleFS;


% ------------------------------------------------------
hcb_3.Label.Position(1) = 3.5;

% ------------------------------------------------------
% Add letter labeling
    
%
xrect_1 = [-3.5, -3.5, -2.6, -2.6, -3.5];
yrect_1 = [1480, 1400, 1400, 1480, 1480];
%
fill(haxs_FCTD, xrect_1, yrect_1, 'w')

%
xrect_2 = 42.6 + [0, 0, 1.5/24, 1.5/24, 0];
yrect_2 = [1500, 1400, 1400, 1500, 1500];
%
for i = 1:length(haxs_moor)
    fill(haxs_moor(i), xrect_2, yrect_2, 'w')
end

%
txtFS = 18;

%
text(haxs_FCTD, -3.35, 1435, 'a)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs_moor(1), 42.61, 1445, 'b)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs_moor(2), 42.61, 1445, 'c)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs_moor(3), 42.61, 1445, 'd)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')


% Add omegaM2 rectangle
%
txtFS = 18;

%
square_1.xlims = [-3.5, -2.3];
square_1.ylims = [1530, 1605]+100;
%
fill(haxs_FCTD, [square_1.xlims(1), square_1.xlims(1), square_1.xlims(2), square_1.xlims(2), square_1.xlims(1)], ...
                [square_1.ylims(2), square_1.ylims(1), square_1.ylims(1), square_1.ylims(2), square_1.ylims(2)], 'w'); 
%
text(haxs_FCTD, mean(square_1.xlims)-0.5, mean(square_1.ylims)-10, '$\mathbf{\omega_{\mathrm{M}_2}}$', 'Interpreter', 'Latex', 'FontSize', 20, 'HorizontalAlignment', 'left')


%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure07.pdf'), 'Resolution', 300)


