%% Attempt to plot both the model and NT1 observations
% on the same figure, which is likely a bit tricky, but
% let's see if I can make something nice.


clear
close all


%%
% ---------------------------------------------------------
% ----------------------- LOAD DATA -----------------------
% ---------------------------------------------------------


%%

%
modelsection = load(fullfile(paper_directory(), 'data', 'extradata', 'data_modelsection.mat'));
%
modelsection = JK_slicebotdepth(modelsection);
             

%% Calculate temperature difference (as a proxy for displacement)

%
modelsection.Tmean = squeeze(nanmean(modelsection.T, 1));

%
modelsection.Tanom = modelsection.T;

%
for i = 1:size(modelsection.T, 1)
    
	%
    modelsection.Tanom(i, :, :) = squeeze(modelsection.Tanom(i, :, :)) - ...
                               modelsection.Tmean;
end


%%

modelsection.T(modelsection.Tmean==0) = NaN;
modelsection.Tanom(modelsection.Tmean==0) = NaN;
modelsection.Tmean(modelsection.Tmean==0) = NaN;


%% NT1 data

%
NT1data = load(fullfile(paper_directory(), 'data', 'shipboard', 'NT1', ...
                                           'NT1_towyodata.mat'));
NT1data = NT1data.NT1data;

%%
% ------------------------------------------------------------
% ---------------- GET INFO ON PROFILE SPACING ---------------
% ------------------------------------------------------------

%%

%
ind_repeat = 4;
ind_z = 300;

%
points.lon = NT1data.CTDdata.transects(ind_repeat).lon(ind_z, :);
points.lat = NT1data.CTDdata.transects(ind_repeat).lat(ind_z, :);
points.time = NT1data.CTDdata.transects(ind_repeat).time(ind_z, :);
%
points.disp = gsw_distance(points.lon, points.lat);
%
points.shipvel = (points.disp ./ diff(points.time))./(24*3600);


%%
% ------------------------------------------------------------
% ------------------------ TRACE RAYS ------------------------
% ------------------------------------------------------------

%% Load stratification

%
modelN2 = load(fullfile(paper_directory(), 'data', 'extradata', 'data_modelstratification.mat'));
modelN2 = modelN2.ray;

% Smooth climatological stratification
%
modelN2.N2_smooth_1 = obmBinAvg(modelN2.z, modelN2.N2, 200, [], @rectwin);
modelN2.N2_smooth_2 = obmBinAvg(modelN2.z, modelN2.N2_smooth_1, 200, [], @rectwin);

% Assign to other variable
profN2.z = modelN2.z;
profN2.dz = profN2.z(2) - profN2.z(1);
profN2.N2 = modelN2.N2_smooth_2;

%
profN2.N2 = profN2.N2(:);


%%

%
f0 = gsw_f(-41.3349);
%
M2_freqrps = 2*pi*(1/12.4206)*(1/3600);

%
M2_freq_aux = M2_freqrps;
M4_freq_aux = 2*M2_freq_aux;
M6_freq_aux = 3*M2_freq_aux;
%
iwM2Char = sqrt((M2_freq_aux^2 - f0^2) ./ (profN2.N2 - M2_freq_aux^2));
iwM4Char = sqrt((M4_freq_aux^2 - f0^2) ./ (profN2.N2 - M4_freq_aux^2));
%
iwM6Char = sqrt((M6_freq_aux^2 - f0^2) ./ (profN2.N2 - M6_freq_aux^2));
iwM6Char = real(iwM6Char); 
            

%% Calculate more specific rays for the plot I want to make
                    
% -----------------------------------
%
plotRays.wvfreq = M2_freq_aux;

% -----------------------------------
%
plotRays.xy0 = [16302.5, 1789.2];

% -----------------------------------
% Leftward (-) or rightward propagation (+)
plotRays.signRLprop = -1;

% -----------------------------------
% Upward (-) or downward propagation (+)
plotRays.signUDprop = -1;

% -----------------------------------
plotRays.dxstep = 10;
%
plotRays.nsteps = 1200;


% -----------------------------------
%
plotRays.iwcharprof = iwM2Char;
%
plotRays.xpts(1) = plotRays.xy0(1);
plotRays.ypts(1) = plotRays.xy0(2);

%
nsteps_aux = plotRays.nsteps;

%
for i2 = 1:nsteps_aux

    %
    angTan_aux = interp1(profN2.z, plotRays.iwcharprof, plotRays.ypts(i2));

    %
    plotRays.xpts(i2+1) = plotRays.xpts(i2) + (plotRays.signRLprop * plotRays.dxstep);

    %
    plotRays.ypts(i2+1) = plotRays.ypts(i2) + (plotRays.signUDprop * (plotRays.dxstep * angTan_aux));
end

%
plotRays.xpts = plotRays.xpts ./ 1000;

% -----------------------------------
%
plotRays.indhaxsplt = [1, 3];



%%
% -------------------------------------------------------
% -------------------------------------------------------
% ---- A VERY SIMPLE PLOT TO CHECK IF IT IS WORTH IT ----
% -------------------------------------------------------
% -------------------------------------------------------

%% Compute x/y grid for NT1 data and put
% these and the model on the same grid

%
lonlat0 = [149.002, -41.3247];

%
for i = 1:NT1data.nrepeats
    %
    [NT1data.CTDdata.transects(i).x, ...
     NT1data.CTDdata.transects(i).y] = ...
                        lonlat2kmgrid(lonlat0, ...
                                      NT1data.CTDdata.transects(i).lon, ...
                                      NT1data.CTDdata.transects(i).lat);
                                  
    %
    [NT1data.LADCPdata.transects(i).x, ...
     NT1data.LADCPdata.transects(i).y] = ...
                        lonlat2kmgrid(lonlat0, ...
                                      NT1data.LADCPdata.transects(i).lon, ...
                                      NT1data.LADCPdata.transects(i).lat);
end

%
[NT1data.bathy.x, NT1data.bathy.y] = ...
                    lonlat2kmgrid(lonlat0, NT1data.bathy.lon, ...
                                           NT1data.bathy.lat);

                              


%%
% -----------------------------------------------------------
% -----------------------------------------------------------
% ----------------------- MAKE FIGURE -----------------------
% -----------------------------------------------------------
% -----------------------------------------------------------
    
%
inds_JK_plt = [6, 9, 12, 2];
inds_NT1_plt = 6:9;    % these kind of must be sequential

%
x0centerplt = 17;

%
Tctrs_aux = 19:0.15:24;
T_ctr_bold = Tctrs_aux(7) .* [1, 1];

%
Tctrs_NT1_aux = 2.2:0.05:3.2;
% % T_ctr_bold = Tctrs_aux(7) .* [1, 1];

%
clrmap.vel = load(fullfile(paper_directory(), 'figures', 'colormaps', 'colormap_velocity.mat'));
clrmap.vel = clrmap.vel.cmapredblue;
    

%%

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.1477    0.1458    0.5650    0.4433];
%
haxs = makeSubPlots(0.1, 0.1, 0.01, ...
                    0.125, 0.125, 0.025, 4, 2);
hold(haxs, 'on')
    
%
for i = 1:length(inds_JK_plt)    % must be 4

    % -----------------------------------
    % Plot JK model

    %
    pcolor(haxs(i), modelsection.x - x0centerplt, modelsection.z, ...
                    squeeze(modelsection.U(inds_JK_plt(i), :, :)));                        
    %
    contour(haxs(i), modelsection.x - x0centerplt, modelsection.z, ...
                     squeeze(modelsection.T(inds_JK_plt(i), :, :)), ...
                     Tctrs_aux, 'k')                         
    %
    contour(haxs(i), modelsection.x - x0centerplt, modelsection.z, ...
                     squeeze(modelsection.T(inds_JK_plt(i), :, :)), ...
                     T_ctr_bold, 'k', 'LineWidth', 2)

    % -----------------------------------
    % Plot NT1 observations
    %
    pcolor(haxs(i+4), NT1data.LADCPdata.transects(inds_NT1_plt(i)).x, ...
                      -NT1data.LADCPdata.z, ...
                      NT1data.LADCPdata.transects(inds_NT1_plt(i)).u)

    %
    contour(haxs(i+4), NT1data.CTDdata.transects(inds_NT1_plt(i)).x, ...
                       -repmat(NT1data.CTDdata.z, 1, NT1data.CTDdata.transects(inds_NT1_plt(i)).nprofiles), ...
                       NT1data.CTDdata.transects(inds_NT1_plt(i)).theta, ...
                       Tctrs_NT1_aux, 'k') 


end
 
        
%
for i = 1:length(haxs)
    %
    shading(haxs(i), 'flat')
    %
    colormap(haxs(i), clrmap.vel)
end



%
xcb = 0.905;
wcb = 0.012;
lblclbFS = 18;
%
hcb_1 = colorbar(haxs(4));
    hcb_1.FontSize = 12;
    hcb_1.Position(1) = xcb;
    hcb_1.Position(3) = wcb;
    hcb_1.Label.String = '$u$ [m s$^{-1}$]';
    hcb_1.Label.Interpreter = 'Latex';
    hcb_1.Label.FontSize = lblclbFS;
%
hcb_2 = colorbar(haxs(8));
    hcb_2.FontSize = 12;
    hcb_2.Position(1) = xcb;
    hcb_2.Position(3) = wcb;
    hcb_2.Label.String = '$u$ [m s$^{-1}$]';
    hcb_2.Label.Interpreter = 'Latex';
    hcb_2.Label.FontSize = lblclbFS;
    hcb_2.Ticks = -0.3:0.15:0.3;


% ------------------------------------------------
% Plot ray
clr_aux = [0, 0, 0];

%
for i = 1:4
    %    
    plot(haxs(i), plotRays(1).xpts - x0centerplt, ...
                  -plotRays(1).ypts, ...
                  '--', 'Color', clr_aux, ...
                  'LineWidth', 4)
end
              
    
% ------------------------------------------------
%
set(haxs, 'FontSize', 12, 'Box', 'on', ...
          'YLim', [-2150, -1380], ...
          'YTick', -2100:200:-1300, ...
          'YTickLabel', {'2100', '1900', '1700', '1500', '1300'})
set(haxs, 'Color', 0.8.*[1, 1, 1])

%
set(haxs(1:4), 'CLim', 0.4.*[-1, 1])
set(haxs(5:8), 'CLim', 0.3.*[-1, 1])
%
set(haxs(1:4), 'XLim', [-5, 3.5])
set(haxs(5:8), 'XLim', [-5, 3.5])
%
set(haxs(1:4), 'XTickLabel', [])
set(haxs([2:4, 6:8]), 'YTickLabel', [])


% ------------------------------------------------

%
for i = 5:length(haxs)
    %
    overlayline(haxs(i), 'v', [15.3679, 17.8271] - x0centerplt, ...
                         '-k', 'LineWidth', 1)
end


% ------------------------------------------------
% Add bathymetry

%
for i = 1:4

    %
    axes(haxs(i))
        fill2Dbathy(modelsection.x - x0centerplt, modelsection.bottomdepth, -5000);

    %
    axes(haxs(i+4))
        %
        fill2Dbathy(NT1data.bathy.x, ...
                    -NT1data.bathy.depth, -5000);
end

% ------------------------------------------------

%
lblFS = 16;
%
ylabel(haxs(1), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(5), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
%
for i = 5:8
    xlabel(haxs(i), '[km]', 'Interpreter', 'Latex', 'FontSize', lblFS)
end

%
title(haxs(1), 'Onshore-offshore reversal', 'Interpreter', 'Latex', 'FontSize', 16)
title(haxs(2), 'Peak offshore flow', 'Interpreter', 'Latex', 'FontSize', 16)
title(haxs(3), 'Offshore-onshore reversal', 'Interpreter', 'Latex', 'FontSize', 16)
title(haxs(4), 'Peak onshore flow', 'Interpreter', 'Latex', 'FontSize', 16)

% ------------------------------------------------------------

%
axes(haxs(2))
    
%
square_1.xlims = [12, 14] - x0centerplt;
square_1.ylims = -[1710, 1820];
%
square_2.xlims = [11.5, 13.8] - x0centerplt + 4.5;
square_2.ylims = -[1380, 1500];

    %
    hf_aux = fill([square_1.xlims(1), square_1.xlims(1), square_1.xlims(2), square_1.xlims(2), square_1.xlims(1)], ...
                  [square_1.ylims(2), square_1.ylims(1), square_1.ylims(1), square_1.ylims(2), square_1.ylims(2)], 'w');
        hf_aux.FaceColor = 1.*[1, 1, 1];
        
	%
    htxt_aux = text(mean(square_1.xlims), mean(square_1.ylims)+15, ...
                                          '$\mathbf{\omega_{\mathrm{M}_2}}$');
        htxt_aux.HorizontalAlignment = 'center';
        htxt_aux.Interpreter = 'Latex';
        htxt_aux.FontSize = 22;
        htxt_aux.Color = clr_aux;

% ------------------------------------------------
% Letter and data labels

% Blank square for letter label on upper left
xsqr_let_lbl = [-5, -5, -3.95, -3.95, -5];
ysqr_let_lbl = [-1500, -1380, -1380, -1500, -1500];
%
for i = 1:length(haxs)
    %
    axes(haxs(i))
        hfillletsqr_aux = fill(xsqr_let_lbl, ysqr_let_lbl, 'w');
% %                 hfillletsqr_aux.LineStyle = 'none';
end

% %     %
% %     xletlbl = -4.8;
% %     yletlbl = -1320;
%
xletlbl = -4.9;
yletlbl = -1440;
%
strletlbls = 'a':'h';
%
for i = 1:length(haxs)
    %
    axes(haxs(i))
        %
        htxt_aux = text(xletlbl, yletlbl, [strletlbls(i) ')']);
            %
            htxt_aux.FontSize = 16;
%                 htxt_aux.FontWeight = 'bold';
end

%
%     xsqrpts = [-5, -5, -2, -2, -5];
ysqrpts = [-2100, -2000, -2000, -2100, -2100];
%
for i = 1:length(haxs)
    %
    if i<=4
        xsqrpts = [-5, -5, -2.75, -2.75, -5];
    else
        xsqrpts = [-5, -5, -0.5, -0.5, -5];
    end

    %
    axes(haxs(i))
        hfill_aux = fill(xsqrpts, ysqrpts, 'w');
            hfill_aux.LineStyle = 'none';
end

%
for i = 1:length(haxs)

    %
    if i<=4
        lbl_data_aux = 'Model';
    else
        lbl_data_aux = 'Towyo (NT1)';
    end

    %
    axes(haxs(i))
        txtData_aux = text(-4.9, -2060, lbl_data_aux);
            txtData_aux.FontSize = 15;
            txtData_aux.Interpreter = 'Latex';

end


% ------------------------------------------------------------
%
moorlbl.T1.x0 = -2.5;
moorlbl.T1.y0 = -1450;
moorlbl.T2.x0 = 1.7;
moorlbl.T2.y0 = -1450;
moorlbl.FS = 18;
%
xlensqr = 1.2;
ylensqr = 120;
%
for i = {'T1', 'T2'}
    %
    moorlbl.(i{1}).xsqr = [(moorlbl.(i{1}).x0 - xlensqr/2), (moorlbl.(i{1}).x0 - xlensqr/2), (moorlbl.(i{1}).x0 + xlensqr/2), (moorlbl.(i{1}).x0 + xlensqr/2), (moorlbl.(i{1}).x0 - xlensqr/2)];
    moorlbl.(i{1}).ysqr = [(moorlbl.(i{1}).y0 - ylensqr/2), (moorlbl.(i{1}).y0 + ylensqr/2), (moorlbl.(i{1}).y0 + ylensqr/2), (moorlbl.(i{1}).y0 - ylensqr/2), (moorlbl.(i{1}).y0 - ylensqr/2)];
end
%
axes(haxs(5))
    %
    hfsqr1 = fill(moorlbl.T1.xsqr, moorlbl.T1.ysqr, 'w');
    hfsqr2 = fill(moorlbl.T2.xsqr, moorlbl.T2.ysqr, 'w');

    %
    htxt_T1 = text(moorlbl.T1.x0, moorlbl.T1.y0 - 10, 'T1');
    htxt_T2 = text(moorlbl.T2.x0, moorlbl.T2.y0 - 10, 'T2');
    %
    htxt_T1.FontSize = moorlbl.FS;
    htxt_T2.FontSize = moorlbl.FS;
    %
    htxt_T1.HorizontalAlignment = 'center';
    htxt_T2.HorizontalAlignment = 'center';


%
set(haxs, 'YLim', [-2100, -1380])

%
hcb_1.Label.Position(1) = hcb_1.Label.Position(1) + 0.5;    % to align with the label on the other colorbar
    

%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure05.pdf'), 'Resolution', 300)



