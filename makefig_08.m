%% Plot tidal ellipses with dissipation

clear
close all


%% Load data

%
dataload.NT1 = load(fullfile(paper_directory(), 'data', 'extradata', 'data_ellipse_NT1.mat'));
dataload.T1 = load(fullfile(paper_directory(), 'data', 'extradata', 'data_ellipse_T1.mat'));
dataload.T2 = load(fullfile(paper_directory(), 'data', 'extradata', 'data_ellipse_T2.mat'));

%
dataload.NT1 = dataload.NT1.NT1fit;
dataload.T1 = dataload.T1.ellipsedata;
dataload.T2 = dataload.T2.ellipsedata;


%% Mooring locations 

mooringlocs = load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));
mooringlocs = mooringlocs.mooringlocs;


%% Load bathymetry

%
bathyTasm = loadbathyTTIDE([148.8, 149.2, -41.4, -41.25]);

%
bathyTasm.depth = -double(bathyTasm.depth);

%
[bathyTasm.longrid, ...
 bathyTasm.latgrid] = meshgrid(bathyTasm.lon, bathyTasm.lat);

%
lonlat0 = dataload.NT1.lonlat0;

%
[bathyTasm.xgrid, ...
 bathyTasm.ygrid] = lonlat2kmgrid(lonlat0, ...
                                  bathyTasm.longrid, bathyTasm.latgrid);
                              

%%
% ----------------------------------
% -- DEFINE VARIABLES FOR PLOTTING -
% ----------------------------------


%%

M2period = 1./dataload.NT1.M2freq;


%% T1

%
nrows = dataload.T1.phasegridded.Nptspercycle;
ncols = length(dataload.T1.phasegridded.u)/nrows;
%
u_array = reshape(dataload.T1.phasegridded.u, nrows, ncols);
v_array = reshape(dataload.T1.phasegridded.v, nrows, ncols);
d_array = reshape(dataload.T1.phasegridded.epsilon, nrows, ncols);
%
u_T1 = mean(u_array(:, 1:6), 2);
v_T1 = mean(v_array(:, 1:6), 2);
d_T1 = mean(d_array(:, 1:6), 2);

% Relative distances (integrate)
dt_int = (24*3600*M2period)/nrows;
xdisp_T1 = [0; cumsum(u_T1(1:end-1)).*dt_int];
ydisp_T1 = [0; cumsum(v_T1(1:end-1)).*dt_int];


%% T2

%
nrows = dataload.T2.phasegridded.Nptspercycle;
ncols = length(dataload.T2.phasegridded.u)/nrows;
%
u_array = reshape(dataload.T2.phasegridded.u, nrows, ncols);
v_array = reshape(dataload.T2.phasegridded.v, nrows, ncols);
d_array = reshape(dataload.T2.phasegridded.epsilon, nrows, ncols);
%
u_T2 = mean(u_array, 2);
v_T2 = mean(v_array, 2);
d_T2 = mean(d_array, 2);

% Relative distances (integrate)
dt_int = (24*3600*M2period)/nrows;
xdisp_T2 = [0; cumsum(u_T2(1:end-1)).*dt_int];
ydisp_T2 = [0; cumsum(v_T2(1:end-1)).*dt_int];


%% Do harmonic fit for moored dissipation

% ---------------------------------------------------------
% T1
%
npts_T1 = length(d_T1);
%
phase_aux = (0 : (npts_T1-1)).*(2*pi./npts_T1);
phase_aux = phase_aux(:);
%
yepsi_aux = log10(d_T1(:));

%
x_array = [ones(length(phase_aux), 1), phase_aux, ...
           cos(phase_aux), ...
           sin(phase_aux)];
    
%
mcoefs_epsilon_T1 = regress(yepsi_aux, x_array);

%
d_T1_fit = 10.^(mcoefs_epsilon_T1(1) + ...
                mcoefs_epsilon_T1(3) .* cos(phase_aux) + ...
                mcoefs_epsilon_T1(4) .* sin(phase_aux));
%
d_T1 = d_T1_fit;


% ---------------------------------------------------------
% T2
%
npts_T2 = length(d_T2);
%
phase_aux = (0 : (npts_T2-1)).*(2*pi./npts_T2);
phase_aux = phase_aux(:);
%
yepsi_aux = log10(d_T2(:));

%
x_array = [ones(length(phase_aux), 1), phase_aux, ...
           cos(phase_aux), ...
           sin(phase_aux)];
    
%
mcoefs_epsilon_T2 = regress(yepsi_aux, x_array);

%
d_T2_fit = 10.^(mcoefs_epsilon_T2(1) + ...
                mcoefs_epsilon_T2(3) .* cos(phase_aux) + ...
                mcoefs_epsilon_T2(4) .* sin(phase_aux));
%
d_T2 = d_T2_fit;


%% NT1 -- onshore ellipse

%
u_NT1_left = dataload.NT1.section(1).ellipse.u(:);
v_NT1_left = dataload.NT1.section(1).ellipse.v(:);
d_NT1_left = dataload.NT1.section(1).ellipse.epsilon(:);

%
dt_int = (24*3600*M2period)/length(dataload.NT1.section(1).ellipse.u);
xdisp_NT1_left = [0; cumsum(u_NT1_left(1:end-1)).*dt_int];
ydisp_NT1_left = [0; cumsum(v_NT1_left(1:end-1)).*dt_int];


%% NT1 -- offshore ellipse

%
u_NT1_right = dataload.NT1.section(2).ellipse.u(:);
v_NT1_right = dataload.NT1.section(2).ellipse.v(:);
d_NT1_right = dataload.NT1.section(2).ellipse.epsilon(:);

%
dt_int = (24*3600*M2period)/length(dataload.NT1.section(2).ellipse.u);
xdisp_NT1_right = [0; cumsum(u_NT1_right(1:end-1)).*dt_int];
ydisp_NT1_right = [0; cumsum(v_NT1_right(1:end-1)).*dt_int];


%% Repeat value at 0, 2*pi so ellipses are closed

%
xdisp_T1 = [xdisp_T1; xdisp_T1(1)];
ydisp_T1 = [ydisp_T1; ydisp_T1(1)];
d_T1 = [d_T1; d_T1(1)];

%
xdisp_T2 = [xdisp_T2; xdisp_T2(1)];
ydisp_T2 = [ydisp_T2; ydisp_T2(1)];
d_T2 = [d_T2; d_T2(1)];

%
xdisp_NT1_left = [xdisp_NT1_left; xdisp_NT1_left(1)];
ydisp_NT1_left = [ydisp_NT1_left; ydisp_NT1_left(1)];
d_NT1_left = [d_NT1_left; d_NT1_left(1)];

%
xdisp_NT1_right = [xdisp_NT1_right; xdisp_NT1_right(1)];
ydisp_NT1_right = [ydisp_NT1_right; ydisp_NT1_right(1)];
d_NT1_right = [d_NT1_right; d_NT1_right(1)];


%% Get reference (x, y) to georefrence data

%
[xloc_T1, yloc_T1] = lonlat2kmgrid(lonlat0, ...
                                   mooringlocs.T1.longitude, ...
                                   mooringlocs.T1.latitude);
%
[xloc_T2, yloc_T2] = lonlat2kmgrid(lonlat0, ...
                                   mooringlocs.T2.longitude, ...
                                   mooringlocs.T2.latitude);
%
xloc_NT1_left = mean(dataload.NT1.section(1).x);
yloc_NT1_left = mean(dataload.NT1.section(1).y);
%
xloc_NT1_right = mean(dataload.NT1.section(2).x);
yloc_NT1_right = mean(dataload.NT1.section(2).y);


%%

%
pltellip.T1.x0 = xloc_T1;
pltellip.T1.y0 = yloc_T1;

%
pltellip.T2.x0 = xloc_T2;
pltellip.T2.y0 = yloc_T2;

%
pltellip.NT1left.x0 = xloc_NT1_left;
pltellip.NT1left.y0 = yloc_NT1_left;

%
pltellip.NT1right.x0 = xloc_NT1_right;
pltellip.NT1right.y0 = yloc_NT1_right;



%% Center ellipses in real (x, y) coordinates

% Define scales for the mean flow vectors, so
% that the ellipses are scales consistently
%
scale_mag = 0.1;
scale_arrow = 12;
%
scale_factor = scale_mag * scale_arrow;


%% Scale as above and get mean position

%
xdisp_T1 = scale_factor*xdisp_T1;
ydisp_T1 = scale_factor*ydisp_T1;
%
xdisp_T2 = scale_factor*xdisp_T2;
ydisp_T2 = scale_factor*ydisp_T2;
%
xdisp_NT1_left = scale_factor*xdisp_NT1_left;
ydisp_NT1_left = scale_factor*ydisp_NT1_left;
%
xdisp_NT1_right = scale_factor*xdisp_NT1_right;
ydisp_NT1_right = scale_factor*ydisp_NT1_right;

%
xavg_disp_T1 = mean(xdisp_T1);    yavg_disp_T1 = mean(ydisp_T1);
xavg_disp_T2 = mean(xdisp_T2);    yavg_disp_T2 = mean(ydisp_T2);
%
xavg_disp_NT1_left = mean(xdisp_NT1_left);    yavg_disp_NT1_left = mean(ydisp_NT1_left);
xavg_disp_NT1_right = mean(xdisp_NT1_right);    yavg_disp_NT1_right = mean(ydisp_NT1_right);


%%

%
pltellip.T1.x = xdisp_T1 - xavg_disp_T1 + 1000*xloc_T1;
pltellip.T1.y = ydisp_T1 - yavg_disp_T1 + 1000*yloc_T1;
%
pltellip.T1.epsilon = d_T1;

%
pltellip.T2.x = xdisp_T2 - xavg_disp_T2 + 1000*xloc_T2;
pltellip.T2.y = ydisp_T2 - yavg_disp_T2 + 1000*yloc_T2;
%
pltellip.T2.epsilon = d_T2;

%
pltellip.NT1left.x = xdisp_NT1_left - xavg_disp_NT1_left + 1000*xloc_NT1_left;
pltellip.NT1left.y = ydisp_NT1_left - yavg_disp_NT1_left + 1000*yloc_NT1_left;
%
pltellip.NT1left.epsilon = d_NT1_left;

%
pltellip.NT1right.x = xdisp_NT1_right - xavg_disp_NT1_right + 1000*xloc_NT1_right;
pltellip.NT1right.y = ydisp_NT1_right - yavg_disp_NT1_right + 1000*yloc_NT1_right;
%
pltellip.NT1right.epsilon = d_NT1_right;

                    
%% Quick check

% % %
% % figure
% %     plot(xdisp_T1, ydisp_T1, '.-k')
% %     hold on
% %     plot(xdisp_T2, ydisp_T2, '.-b')
% %     
% %     %
% %     plot(xdisp_NT1_left, ydisp_NT1_left, '.-m')
% %     plot(xdisp_NT1_right, ydisp_NT1_right, '.-r')
% %     
% %     axis equal
    
    
% % %
% % figure
% %     %
% %     hold on
% %     %
% %     plot(xdisp_T1, ydisp_T1, '-k', 'LineWidth', 6)
% %     plot(xdisp_T2, ydisp_T2, '-k', 'LineWidth', 6)
% %     %
% %     plot(xdisp_NT1_left, ydisp_NT1_left, '-k', 'LineWidth', 6)
% %     plot(xdisp_NT1_right, ydisp_NT1_right, '-k', 'LineWidth', 6)
% %     
% %     %
% %     mkSZ = 62;
% %     %
% %     scatter(xdisp_T1, ydisp_T1, mkSZ, log10(d_T1), 'filled')
% %     scatter(xdisp_T2, ydisp_T2, mkSZ, log10(d_T2), 'filled')
% %     %
% %     scatter(xdisp_NT1_left, ydisp_NT1_left, mkSZ, log10(d_NT1_left), 'filled')
% %     scatter(xdisp_NT1_right, ydisp_NT1_right, mkSZ, log10(d_NT1_right), 'filled')
% %     
% %     %
% %     colorbar
% %     
% %     %
% %     axis equal
    
    
% % %
% % figure
% %     %
% %     hold on
% %     %
% %     plot(pltellip.T1.x, pltellip.T1.y, '-k', 'LineWidth', 6)
% %     plot(pltellip.T2.x, pltellip.T2.y, '-k', 'LineWidth', 6)
% %     %
% %     plot(pltellip.NT1left.x, pltellip.NT1left.y, '-k', 'LineWidth', 6)
% %     plot(pltellip.NT1right.x, pltellip.NT1right.y, '-k', 'LineWidth', 6)
% %     
% %     %
% %     mkSZ = 62;
% %     %
% %     scatter(pltellip.T1.x, pltellip.T1.y, mkSZ, log10(d_T1), 'filled')
% %     scatter(pltellip.T2.x, pltellip.T2.y, mkSZ, log10(d_T2), 'filled')
% %     %
% %     scatter(pltellip.NT1left.x, pltellip.NT1left.y, mkSZ, log10(d_NT1_left), 'filled')
% %     scatter(pltellip.NT1right.x, pltellip.NT1right.y, mkSZ, log10(d_NT1_right), 'filled')
% %     
% %     %
% %     colorbar
% %     
% %     %
% %     axis equal
    
      
%% Get mean flow

%
pltellip.T1.meanflow = dataload.T1.meanflow;
pltellip.T2.meanflow = dataload.T2.meanflow;
%
pltellip.NT1left.meanflow = dataload.NT1.section(1).meanflow;
pltellip.NT1right.meanflow = dataload.NT1.section(2).meanflow;


%%
% -----------------------------------------------------------
% -----------------------------------------------------------
% ----------------------- MAKE FIGURE -----------------------
% -----------------------------------------------------------
% -----------------------------------------------------------


%%

%
xaxs_lims = [-2.4, 2];
yaxs_lims = [-1.5, 3.25];

%
cmapbathy = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_bathymetry.mat'));
cmapbathy = cmapbathy.cmapblues;
%
cmapepsi = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_epsilon.mat'));
cmapepsi = cmapepsi.cmaporange;


%% Plot figure

%
hfig = figure;
hfig.Units = 'normalized';
%
% % hfig.Position
%
%
haxs_backg = makeSubPlots(0.1, 0.3, 0.1, ...
                          0.1, 0.2, 0.1, 1, 1);
hold(haxs_backg, 'on')

% ------------------------------------
    %
    hpc = pcolor(haxs_backg, bathyTasm.xgrid, bathyTasm.ygrid, bathyTasm.depth);
	%
    %
    contour(haxs_backg, bathyTasm.xgrid, bathyTasm.ygrid, bathyTasm.depth, 1500:20:2400, 'k')
    %
    contour(haxs_backg, bathyTasm.xgrid, bathyTasm.ygrid, bathyTasm.depth, 1000:100:4000, 'k', 'LineWidth', 4)

    %
    plot(haxs_backg, dataload.NT1.section(1).x(:), dataload.NT1.section(1).y(:), '.k')
    plot(haxs_backg, dataload.NT1.section(2).x(:), dataload.NT1.section(2).y(:), '.k')
    
    %
    lnWD_aux = 14;
    %
    plot(haxs_backg, pltellip.T1.x./1000, pltellip.T1.y./1000, 'k', 'LineWidth', lnWD_aux)
    plot(haxs_backg, pltellip.T2.x./1000, pltellip.T2.y./1000, 'k', 'LineWidth', lnWD_aux)
    %
    plot(haxs_backg, pltellip.NT1left.x./1000, pltellip.NT1left.y./1000, 'k', 'LineWidth', lnWD_aux)
    plot(haxs_backg, pltellip.NT1right.x./1000, pltellip.NT1right.y./1000, 'k', 'LineWidth', lnWD_aux)
    
    % Plot same as above, but changing the order (break point), so that
    % Matlab plot continuous lines
    %
    plot(haxs_backg, pltellip.T1.x([10:25, 26:33, 1:9])./1000, pltellip.T1.y([10:25, 26:33, 1:9])./1000, 'k', 'LineWidth', lnWD_aux)
    plot(haxs_backg, pltellip.T2.x([10:25, 26:33, 1:9])./1000, pltellip.T2.y([10:25, 26:33, 1:9])./1000, 'k', 'LineWidth', lnWD_aux)
    %
    plot(haxs_backg, pltellip.NT1left.x([10:25, 26:33, 1:9])./1000, pltellip.NT1left.y([10:25, 26:33, 1:9])./1000, 'k', 'LineWidth', lnWD_aux)
    plot(haxs_backg, pltellip.NT1right.x([10:25, 26:33, 1:9])./1000, pltellip.NT1right.y([10:25, 26:33, 1:9])./1000, 'k', 'LineWidth', lnWD_aux)

    % Markers for T1 and T2 mooring
    mkszmoor = 90;
    plot(pltellip.T1.x0, pltellip.T1.y0, '.k', 'MarkerSize', mkszmoor);
    plot(pltellip.T2.x0, pltellip.T2.y0, '.k', 'MarkerSize', mkszmoor);
    
% ------------------------------------
%
haxs_epsi = axes('Position', haxs_backg.Position);
hold(haxs_epsi, 'on')

    
    %
    epsi_scatter_SZ = 100;
    
    %
    scatter(haxs_epsi, pltellip.T1.x./1000, pltellip.T1.y./1000, epsi_scatter_SZ, log10(pltellip.T1.epsilon), 'filled')
    scatter(haxs_epsi, pltellip.T2.x./1000, pltellip.T2.y./1000, epsi_scatter_SZ, log10(pltellip.T2.epsilon), 'filled')
                                  
    %
    scatter(haxs_epsi, pltellip.NT1left.x./1000, pltellip.NT1left.y./1000, epsi_scatter_SZ, log10(pltellip.NT1left.epsilon), 'filled')
    scatter(haxs_epsi, pltellip.NT1right.x./1000, pltellip.NT1right.y./1000, epsi_scatter_SZ, log10(pltellip.NT1right.epsilon), 'filled')
           
    


% ------------------------------------

%
ycb = 0.2;
xwcb = 0.02;
ywcb = 0.7;
%
hcb = colorbar(haxs_epsi);
    hcb.Position = [0.66, ycb, xwcb, ywcb];
    hcb.Label.Interpreter = 'Latex';
    hcb.Label.String = '$\log_{10}\epsilon~[\mathrm{m}^2~\mathrm{s}^{-3}]$';
    hcb.FontSize = 16;
% %     hcb.Ticks = -10 : 1 : -7;
    hcb.Ticks = -8.5 : 0.5 : -7.5;
    %
    hcb.Label.FontSize = 18;

    
        
% ------------------------------------

    
%
hpc.LineStyle = 'none';
hpc.FaceAlpha = 0.6;

%
colormap(haxs_backg, cmapbathy);
colormap(haxs_epsi, cmapepsi);

%
set(haxs_backg, 'DataAspectRatio', [1, 1, 1], 'XLim', xaxs_lims, 'YLim', yaxs_lims)
set(haxs_backg, 'CLim', [1750, 2100])

%
set(haxs_epsi, 'DataAspectRatio', [1, 1, 1], 'XLim', xaxs_lims, 'YLim', yaxs_lims)
set(haxs_epsi, 'CLim', [-8.5, -7.4])
    
%
set(haxs_backg, 'FontSize', 14, 'XTick', -3:1:3, 'YTick', -1:1:4)

%
haxs_epsi.Visible = 'off';

% ------------------------------------

%
xlabel(haxs_backg, 'cross-shore [km]', 'Interpreter', 'Latex', 'FontSize', 16)
ylabel(haxs_backg, 'alongshore [km]', 'Interpreter', 'Latex', 'FontSize', 16)

%
haxs_epsi.Visible = 'off';

% ------------------------------------------------------------
% Plot arrows for the mean flow

%
shaftwidth = 0.1;
arrwhdwidth = 0.1;
arrwhdlen = 0.2;
%
xcoordbase = [0, 0, (1-arrwhdlen), (1-arrwhdlen), 1, (1-arrwhdlen), (1-arrwhdlen), 0];
ycoordbase = [-shaftwidth/2, +shaftwidth/2, +shaftwidth/2, +shaftwidth/2 + arrwhdwidth, 0, (-shaftwidth/2 - arrwhdwidth), -shaftwidth/2, -shaftwidth/2];

%
list_data = {'T1', 'T2', 'NT1left', 'NT1right'};

%
for i = 1:length(list_data)
    
    %
    x_arrow_aux = xcoordbase(:) .* scale_arrow * ...
                                   sqrt(pltellip.(list_data{i}).meanflow.u.^2 + ...
                                        pltellip.(list_data{i}).meanflow.v.^2);
	% Change neck position so it does not scale with magnitude
    x_arrow_aux([3:4, 6:7]) = x_arrow_aux(5) - arrwhdlen;
                                    
	%
    y_arrow_aux = ycoordbase(:);
    
    %
    dir_flow = atan2(pltellip.(list_data{i}).meanflow.v, pltellip.(list_data{i}).meanflow.u);
    
    %
    rot_matrix = [cos(dir_flow), -sin(dir_flow); sin(dir_flow), cos(dir_flow)];
    
    %
    xy_rot_aux = rot_matrix * [x_arrow_aux.'; y_arrow_aux.'];
    
    %
    x_arrow_plt = xy_rot_aux(1, :) + pltellip.(list_data{i}).x0;
    y_arrow_plt = xy_rot_aux(2, :) + pltellip.(list_data{i}).y0;
    
    %
    hfill_aux = fill(haxs_epsi, x_arrow_plt, y_arrow_plt, 'g');
    
end

% ------------------------------------------------------------
% Mean flow scale

%
xwidthsqr = 1.5;
yheightsqr = 0.95;

    %
    fill([(xaxs_lims(2) - xwidthsqr), xaxs_lims(2), xaxs_lims(2), (xaxs_lims(2) - xwidthsqr), (xaxs_lims(2) - xwidthsqr)], ...
         [(yaxs_lims(1)+yheightsqr), (yaxs_lims(1)+yheightsqr), yaxs_lims(1), yaxs_lims(1), yaxs_lims(1)], 'w');

    %
    htxtscale = text((xaxs_lims(2) - 0.74), ...
                     yaxs_lims(1) + 0.75 - 0.64, ...
                     ['$' num2str(scale_mag, '%.1f') '~\mathrm{m}~\mathrm{s}^{-1}$']);
        htxtscale.Interpreter = 'Latex';
        htxtscale.FontSize = 14;
        htxtscale.HorizontalAlignment = 'center';
    %
    htxtscale_lbl = text((xaxs_lims(2) - 0.75), ...
                         yaxs_lims(1) + 0.7, 'mean flow');
        htxtscale_lbl.Interpreter = 'Latex';
        htxtscale_lbl.FontSize = 14;
        htxtscale_lbl.HorizontalAlignment = 'center'; 
        
        
%
x_scale = xcoordbase(:) .* scale_arrow * scale_mag;
% Change neck position so it does not scale with magnitude
x_scale([3:4, 6:7]) = x_scale(5) - arrwhdlen;
%
y_scale = ycoordbase(:);
    %
    hfill_scale = fill(haxs_epsi, x_scale + 0.68, y_scale - 1.075, 'g');

    
% ------------------------------------------------------------
% Annotations for moorings and NT1


% -----------------------
%
txtFS = 16;
%
htxtT1 = text(pltellip.T1.x0, pltellip.T1.y0, 'T1');
    htxtT1.FontSize = txtFS;
    htxtT1.Color = [1, 1, 1];
    htxtT1.HorizontalAlignment = 'center';
    htxtT1.FontWeight = 'bold';
%
htxtT2 = text(pltellip.T2.x0, pltellip.T2.y0, 'T2');
    htxtT2.FontSize = txtFS;
    htxtT2.Color = [1, 1, 1];
    htxtT2.HorizontalAlignment = 'center';
    htxtT2.FontWeight = 'bold';
        
% -----------------------
%
xdbox = 0.88;
ydbox = 0.4;
%
x0NT1 = -0.7;
y0NT1 = 1.6;

%
hfNT1aux = fill([x0NT1, x0NT1, (x0NT1 + xdbox), (x0NT1 + xdbox), xdbox], ...
                [y0NT1, (y0NT1 + ydbox), (y0NT1 + ydbox), y0NT1, y0NT1], 'k');
    hfNT1aux.LineStyle = 'none';

%
htxtNT1_aux = text(mean([x0NT1, (x0NT1 + xdbox)]), ...
                   mean([y0NT1, (y0NT1 + ydbox)]), 'NT1');
    htxtNT1_aux.Color = [1, 1, 1];
    htxtNT1_aux.HorizontalAlignment = 'center';
    htxtNT1_aux.FontWeight = 'bold';
    htxtNT1_aux.FontSize = txtFS;


%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure08.pdf'), 'Resolution', 300)
