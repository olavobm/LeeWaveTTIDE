%% Make T1 timeseries figure

clear
close all


%%
% -------------------------------------------
% ---------------- LOAD DATA ----------------
% -------------------------------------------

% Load velocity
%
load(fullfile(paper_directory(), 'data', 'mooring', 'level_2', 'T1_velocity_L2.mat'))
%
adcpfilt.D2 = load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_velocity_D2.mat'));
adcpfilt.subtidal = load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_velocity_subtidal.mat'));
%
adcpfilt.D2 = adcpfilt.D2.adcpdata;
adcpfilt.subtidal = adcpfilt.subtidal.adcpdata;

% Load turbulence
load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_turbulence.mat'))


%%
% --------------------------------------------
% ------- COMPUTE TIMESERIES FOR PANELS ------
% --------------- (d) and (e) ----------------
% --------------------------------------------

%%

ztopavg = 1600;


%% Depth-average velocity

%
lbelowtop_velD2 = adcpfilt.D2.z >= ztopavg;
lbelowtop_velsubtidal = adcpfilt.subtidal.z >= ztopavg;
% (these two should be equal though)

%
pltvars.velocity.dtime = adcpfilt.D2.dtime;
pltvars.velocity.yday = adcpfilt.D2.yday;
% (should be the same as in subtidal)

%
pltvars.velocity.D2.u = mean(adcpfilt.D2.u(lbelowtop_velD2, :), 1, 'omitnan');
pltvars.velocity.D2.v = mean(adcpfilt.D2.v(lbelowtop_velD2, :), 1, 'omitnan');
%
pltvars.velocity.subtidal.u = mean(adcpfilt.subtidal.u(lbelowtop_velsubtidal, :), 1, 'omitnan');
pltvars.velocity.subtidal.v = mean(adcpfilt.subtidal.v(lbelowtop_velsubtidal, :), 1, 'omitnan');


%% Depth-integrate dissipation

%
lbelowtop_epsi = tchaindata.turbulence.zgrid >= ztopavg;

%
pltvars.turbulence.dtime = tchaindata.turbulence.dtime;
pltvars.turbulence.yday = tchaindata.turbulence.yday;
pltvars.turbulence.epsilon = 1025 * tchaindata.turbulence.dzgrid * ...
                                    sum(tchaindata.turbulence.epsilon_gridded(lbelowtop_epsi, :), 1, 'omitnan');

%
dtime_start = datetime(2015, 01, 18, 14, 00, 00, 'TimeZone', tchaindata.turbulence.dtime.TimeZone);
dtime_stop = datetime(2015, 03, 04, 03, 00, 00, 'TimeZone', tchaindata.turbulence.dtime.TimeZone);
%
dt_step = hours(1);

%
ind_start = find(pltvars.turbulence.dtime >= (dtime_start-(dt_step/2)), 1, 'first');
ind_stop = find(pltvars.turbulence.dtime   < (dtime_stop+(dt_step/2)), 1, 'last');
%
ind_get = ind_start : 1 : ind_stop;

%
Nrows = dt_step/tchaindata.turbulence.dt;
Ncols = length(ind_get) / Nrows;

%
epsilon_array = reshape(pltvars.turbulence.epsilon(ind_get), Nrows, Ncols);


%
pltvars.turbulence.dtime_hourly = dtime_start : dt_step : dtime_stop;
pltvars.turbulence.yday_hourly = datenum(pltvars.turbulence.dtime_hourly) - datenum(2015, 1, 1);
pltvars.turbulence.epsilon_hourly = mean(epsilon_array, 1);

%
pltvars.turbulence.epsilon_12hour = NaN(1, length(pltvars.turbulence.epsilon_hourly));

%
for i = 1:length(pltvars.turbulence.dtime_hourly)
    %
    linwindow = (pltvars.turbulence.dtime >= (pltvars.turbulence.dtime_hourly(i) - hours(6))) & ...
                (pltvars.turbulence.dtime  < (pltvars.turbulence.dtime_hourly(i) + hours(6)));
    %
    pltvars.turbulence.epsilon_12hour(i) = mean(pltvars.turbulence.epsilon(linwindow));
end


%%
% --------------------------------------------
% ---------------- MAKE FIGURE ---------------
% --------------------------------------------


%
xlimall = [17.1736, 62.1875];
yaxslims = [1500, 1920];

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.1455    0.1417    0.4559    0.6250];
%
%
haxs = makeSubPlots(0.1, 0.2, 0.1, ...
                    0.1, 0.125, 0.01, 1, 5);
hold(haxs, 'on')

% --------------------------------------------------------
%
for i = 3:5

    %
    axes(haxs(i))
    %
    xlim(xlimall)
    %
    if i==4
        ylim(0.2.*[-1, 1])
    elseif i==3
        ylim([1350, yaxslims(2)])
    elseif i==5
        ylim([1e-5, 1])
    end

    % 
    time_shading = [18.4, 22; ...
                    42.1978, 46.8];
    %
    clrs_shading = [0, 0.5, 0; ...
                    0, 0, 0.5];
                
    %
    sliceShading(time_shading, clrs_shading, 0.5);
end

% --------------------------------------------------------
%
htxt_1_aux = text(haxs(4), 22.25, -0.105, 'T2 + NT1 towyo');
htxt_2_aux = text(haxs(4), 36, 0.1, 'Fast-CTD');

%
htxt_1_aux.Color = clrs_shading(1, :);
htxt_2_aux.Color = clrs_shading(2, :);

%
FSaux = 16;
%
htxt_1_aux.FontSize = FSaux;
htxt_2_aux.FontSize = FSaux;


% --------------------------------------------------------

    %
    hpc_1 = pcolor(haxs(1), adcpdata.yday, adcpdata.z, adcpdata.u);
    hpc_2 = pcolor(haxs(2), adcpdata.yday, adcpdata.z, adcpdata.v);
    %
    hpc_3 = pcolor(haxs(3), tchaindata.turbulence.yday, ...
                            tchaindata.turbulence.zgrid, ...
                            log10(tchaindata.turbulence.epsilon_gridded));

%
hpc_1.LineStyle = 'none';
hpc_2.LineStyle = 'none';
hpc_3.LineStyle = 'none';


    %
    hp_UD2 = plot(haxs(4), pltvars.velocity.yday, pltvars.velocity.D2.u, 'r', 'LineWidth', 1.6);
    %
    hp_Usub = plot(haxs(4), pltvars.velocity.yday, pltvars.velocity.subtidal.u, 'b', 'LineWidth', 1.5);
    hp_Vsub = plot(haxs(4), pltvars.velocity.yday, pltvars.velocity.subtidal.v, 'k', 'LineWidth', 1.5);

    %
    hp_epsi1 = plot(haxs(5), pltvars.turbulence.yday_hourly, pltvars.turbulence.epsilon_hourly, '.-k');
    hp_epsi2 = plot(haxs(5), pltvars.turbulence.yday_hourly, pltvars.turbulence.epsilon_12hour, '-r', 'LineWidth', 2.5);

%
set(haxs(5)', 'YScale', 'log', 'YLim', [2e-4, 2e-2], 'YTick', [10.^(-5:1:-1)])

    %
    plot(haxs(4), xlimall, [0, 0], '-k')


% --------------------------------------------------
%
cmapdata.vel = load('colormap_velocity.mat');
cmapdata.epsi = load('colormap_epsilon.mat');
%
cmapdata.vel = cmapdata.vel.cmapredblue;
cmapdata.epsi = cmapdata.epsi.cmaporange;

%
set(haxs(1:2), 'CLim', 0.2.*[-1, 1], 'Colormap', cmapdata.vel)
set(haxs(3), 'CLim', [-9, -7], 'Colormap', cmapdata.epsi)

%
xdcb = 0.81;
wcb = 0.015;
%
hcb_1 = colorbar(haxs(2));
    hcb_1.Position(1) = xdcb;
    hcb_1.Position(3) = wcb;
    hcb_1.Position(4) = sum(haxs(1).Position([2, 4])) - haxs(2).Position(2);
    hcb_1.Label.Interpreter = 'Latex';
    hcb_1.Label.FontSize = 18;
    hcb_1.Label.String = '[m s$^{-1}$]';
    hcb_1.Ticks = -0.2 : 0.1 : 0.2;
    hcb_1.FontSize = 14;
%
hcb_2 = colorbar(haxs(3));
    hcb_2.Position(1) = xdcb;
    hcb_2.Position(3) = wcb;
    hcb_2.Label.Interpreter = 'Latex';
    hcb_2.Label.FontSize = 16;
    hcb_2.Label.String = {'$\log_{10}\epsilon$';'$[\mathrm{m}^2\mathrm{s}^{-3}]$'};
    hcb_2.FontSize = 14;

% --------------------------------------------------
%
set(haxs, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'XLim', xlimall)
set(haxs(1:3), 'YLim', [1400, 1920], 'YDir', 'reverse')
set(haxs(4), 'YLim', 0.13.*[-1, 1])
set(haxs(5), 'YLim', [3e-4, 1.5e-2], 'YScale', 'log', 'YTick', [1e-4, 1e-3, 1e-2, 1e-1])
%
for i = 1:length(haxs)
    haxs(i).YAxis.FontSize = 11;
end
%
set(haxs, 'Color', 0.8.*[1, 1, 1])
set(haxs(1:(end-1)), 'XTickLabel', [])

%
xlabel(haxs(end), 'Time [yearday in 2015]', 'Interpreter', 'Latex', 'FontSize', 18)
%
ylblFS = 12;
%
ylabel(haxs(1), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', ylblFS)
ylabel(haxs(2), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', ylblFS)
ylabel(haxs(3), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', ylblFS)
ylabel(haxs(4), '[m s$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', ylblFS)
ylabel(haxs(5), '[W m$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', ylblFS)



% --------------------------------------------------------
% --------------------------------------------------------
    %
% %     haxs_lvel = axes('Position', [0.81, 0.282+0.02, 0.11, 0.147-0.04]);
% %     haxs_lepsi = axes('Position', [0.81, 0.125+0.005, 0.11, 0.147-0.01]);
    haxs_lvel = axes('Position', [0.81, 0.282+0.02, 0.09, 0.147-0.04]);
    haxs_lepsi = axes('Position', [0.81, 0.125+0.005, 0.09, 0.147-0.01]);

        %
        hold(haxs_lvel, 'on');
        hold(haxs_lepsi, 'on');
        %
        haxs_lvel.Box = 'on';
        haxs_lepsi.Box = 'on';

        %
        xaxs_lim_extra = [0, 1];
        yaxs_lim_extra = [0, 1];


        %
        plot(haxs_lvel, [0.03, 0.4], 0.8.*[1, 1], 'Color', hp_UD2.Color, 'LineWidth', 3)
        plot(haxs_lvel, [0.03, 0.4], 0.5.*[1, 1], 'Color', hp_Usub.Color, 'LineWidth', 3)
        plot(haxs_lvel, [0.03, 0.4], 0.2.*[1, 1], 'Color', hp_Vsub.Color, 'LineWidth', 3)

        %
        plot(haxs_lepsi, [0.03, 0.25], 0.75.*[1, 1], 'Color', hp_epsi1.Color, 'LineWidth', 3)
        plot(haxs_lepsi, mean([0.03, 0.25]), 0.75, '.', 'Color', hp_epsi1.Color, 'MarkerSize', 20)
        plot(haxs_lepsi, [0.03, 0.25], 0.25.*[1, 1], 'Color', hp_epsi2.Color, 'LineWidth', 3)

        %
        text(haxs_lvel, 0.55, 0.8, '$U_{\mathrm{D}2}$', 'Interpreter', 'Latex', 'FontSize', 16)
        text(haxs_lvel, 0.55, 0.5, '$\langle U \rangle$', 'Interpreter', 'Latex', 'FontSize', 16)
        text(haxs_lvel, 0.55, 0.2, '$\langle V \rangle$', 'Interpreter', 'Latex', 'FontSize', 16)
% %         %
% %         text(haxs_leg_epsi, 0.4, 0.7, {'hourly';'average'}, 'Interpreter', 'Latex', 'FontSize', 14)
% %         text(haxs_leg_epsi, 0.4, 0.3, {'semidiurnal';'average'}, 'Interpreter', 'Latex', 'FontSize', 14)
        %
        text(haxs_lepsi, 0.3, 0.75+0.08, 'hourly', 'Interpreter', 'Latex', 'FontSize', 14)
        text(haxs_lepsi, 0.3, 0.75-0.08, 'average', 'Interpreter', 'Latex', 'FontSize', 14)
% %         text(haxs_leg_epsi, 0.3, 0.25+0.1, 'semidiurnal', 'Interpreter', 'Latex', 'FontSize', 14)
% %         text(haxs_leg_epsi, 0.3, 0.25-0.1, 'average', 'Interpreter', 'Latex', 'FontSize', 14)
        text(haxs_lepsi, 0.3, 0.25+0.16, 'semi-', 'Interpreter', 'Latex', 'FontSize', 14)
        text(haxs_lepsi, 0.3, 0.25, 'diurnal', 'Interpreter', 'Latex', 'FontSize', 14)
        text(haxs_lepsi, 0.3, 0.25-0.16, 'average', 'Interpreter', 'Latex', 'FontSize', 14)

        %
        haxs_lvel.XLim = xaxs_lim_extra;
        haxs_lvel.YLim = yaxs_lim_extra;
        %
        haxs_lepsi.XLim = xaxs_lim_extra;
        haxs_lepsi.YLim = yaxs_lim_extra;

        %
        haxs_lvel.XTick = [];
        haxs_lvel.YTick = [];
        %
        haxs_lepsi.XTick = [];
        haxs_lepsi.YTick = [];
	% --------------------------------------------
	%
    haxs_period = axes('Position', [haxs(1).Position(1), 0.886, haxs(1).Position(3), 0.08]);
    hold(haxs_period, 'on')
    
        %
        ylim(haxs_period, [0, 1])

        %
        xlimall = xlim(haxs(1));
        %
        pause(1)
        %
        plot(haxs_period, [xlimall(1), 26], 0.3.*[1, 1], '-k')
        plot(haxs_period, [27, 54], 0.3.*[1, 1], '-k')
        plot(haxs_period, [55, xlimall(2)], 0.3.*[1, 1], '-k')
        %
        haxs_period.Visible = 'off';
        
    %
    axes(haxs(1))
        %
        htxt_title = text(27, 1200, 'mooring observations (T1)');
            htxt_title.Interpreter = 'Latex';
            htxt_title.FontSize = 22;
            
        %
        lblperiodFS = 16;
        lblypos = 1225;
        %
        htxt_periodI = text(21.5, 1310, 'period I');
%             htxt_periodI.Interpreter = 'Latex';
            htxt_periodI.FontSize = lblperiodFS;
            htxt_periodI.HorizontalAlignment = 'center';
        %
        htxt_periodII = text(mean(xlimall), 1310, 'period II');
%             htxt_periodII.Interpreter = 'Latex';
            htxt_periodII.FontSize = lblperiodFS;
            htxt_periodII.HorizontalAlignment = 'center';
        %
        htxt_periodIII = text(58.5, 1310, 'period III');
%             htxt_periodIII.Interpreter = 'Latex';
            htxt_periodIII.FontSize = lblperiodFS;
            htxt_periodIII.HorizontalAlignment = 'center';
            
            
    %
    set(haxs_period, 'XLim', xlimall)

    
% --------------------------------------------------------
    
%
xrect_pclr_1 = [17.2, 21, 21, 17.2, 17.2];
%
xrect_epsi_pclr = [17.1736, 20.2, 20.2, 17.1736, 17.1736];
%
xrect_vel = [17.1736, 25, 25, 17.1736, 17.1736];
%
xrect_intepsi = [17.1736, 25, 25, 17.1736, 17.1736];


%
yrect_1_pclr = [1400, 1400, 1530, 1530, 1400];
yrect_2_pclr = [0.08, 0.08, 0.13, 0.13, 0.08];
yrect_3_pclr = [5e-3, 5e-3, 1.5e-2, 1.5e-2, 5e-3];
%

%
fill(haxs(1), xrect_pclr_1, yrect_1_pclr, 'w')
fill(haxs(2), xrect_pclr_1, yrect_1_pclr, 'w')
fill(haxs(3), xrect_pclr_1, yrect_1_pclr, 'w')
%
fill(haxs(4), xrect_vel, yrect_2_pclr, 'w')
fill(haxs(5), xrect_intepsi, yrect_3_pclr, 'w')

%
txtFS = 17;
%
xref = 17.35;
%
text(haxs(1), xref, mean(yrect_1_pclr(2:3)), 'a)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs(1), mean(xrect_pclr_1(1:2))-0.02, mean(yrect_1_pclr(2:3)), '$U$', 'Interpreter', 'Latex', 'FontSize', txtFS)
%
text(haxs(2), xref, mean(yrect_1_pclr(2:3)), 'b)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs(2), mean(xrect_pclr_1(1:2))-0.02, mean(yrect_1_pclr(2:3)), '$V$', 'Interpreter', 'Latex', 'FontSize', txtFS)
%
text(haxs(3), xref, mean(yrect_1_pclr(2:3)), 'c)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs(3), xrect_epsi_pclr(1)+1.75, mean(yrect_1_pclr(2:3)), '$\epsilon$', 'Interpreter', 'Latex', 'FontSize', txtFS+6)
%
text(haxs(4), xref, mean(yrect_2_pclr(2:3)), 'd)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs(4), mean(xrect_vel(1:2))-2.2, mean(yrect_2_pclr(2:3)), 'velocity', 'Interpreter', 'Latex', 'FontSize', txtFS)
%
text(haxs(5), xref, mean(yrect_3_pclr(2:3))-1.5e-3, 'e)', 'FontSize', txtFS, 'HorizontalAlignment', 'left')
text(haxs(5), mean(xrect_intepsi(1:2))-2.25, mean(yrect_3_pclr(2:3))-1.5e-3, '$\int \rho_0 \epsilon~\mathrm{d}z$', 'Interpreter', 'Latex', 'FontSize', txtFS)


%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure04.pdf'), 'Resolution', 300)




