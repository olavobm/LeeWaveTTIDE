%% Make Figure 9 (scatering model).

clear
close all


%%
% --------------------------------------------------------------
% ------------------------- LOAD STUFF -------------------------
% --------------------------------------------------------------

%% Load bathymetry

%
lonbox_slope = [148.632, 150];
latbox_slope = [-41.4, -41.28];
%
indstepload = 1;

%
bathy = loadbathyTTIDE([lonbox_slope, latbox_slope], indstepload);

% Make it positive
bathy.depth = - double(bathy.depth);


%% Mooring location

%
moorloc = load(fullfile(paper_directory(), 'data', 'metadata', 'mooring_locations.mat'));
moorloc = moorloc.mooringlocs;
%
T1loc = moorloc.T1;


%% Load stratification (from Jody's product)

% % %
% % JKN2 = load(['/Users/olavobm/Documents/MATLAB/olavo_research' ...
% %              '/TTIDE/JodysDropbox/TasmaniaRays.mat']);
% % JKN2 = JKN2.ray;

JKN2 = load(fullfile(paper_directory(), 'data', 'extradata', ...
                                        'data_modelstratification.mat'));
JKN2 = JKN2.ray;


% Smooth climatological stratification
%
JKN2.N2_smooth_1 = obmBinAvg(JKN2.z, JKN2.N2, 200, [], @rectwin);
JKN2.N2_smooth_2 = obmBinAvg(JKN2.z, JKN2.N2_smooth_1, 200, [], @rectwin);

% % %
% % newFigDims([5.8889, 7.0972])
% %     plot(JKN2.N2, JKN2.z)
% %     hold on
% % %     plot(JKN2.N2_smooth_1, JKN2.z, 'LineWidth', 2)
% %     plot(JKN2.N2_smooth_2, JKN2.z, 'LineWidth', 2)
% %     axis ij
% %     grid on
% %     set(gca, 'FontSize', 16)
% %     xlabel('N$^2$ [s$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', 22)
% %     ylabel('Depth [m]', 'Interpreter', 'Latex', 'FontSize', 22)
% %     title('Stratification used in CELT', 'Interpreter', 'Latex', 'FontSize', 22)

% Assign to other variable
profN2.z = JKN2.z;
profN2.dz = profN2.z(2) - profN2.z(1);
profN2.N2 = JKN2.N2_smooth_2;



%%
% --------------------------------------------------------------
% ------------------- CALCULATE ON X/Y GRID --------------------
% --------------------------------------------------------------
%
% (make the bump to be x = 0)


%%

%
[bathy.longrid, bathy.latgrid] = meshgrid(bathy.lon, bathy.lat);

%


%%

%
lonlat0 = [149.0024, -41.3394];

%
[bathy.xgrid, bathy.ygrid] = lonlat2kmgrid(lonlat0, ...
                                           bathy.longrid, bathy.latgrid);

                                       
%%
% --------------------------------------------------------------
% -------------------- GET A BOTTOM TRANSECT -------------------
% --------------------------------------------------------------

%%

%
indsecplt = round(length(bathy.lat)/2);

%
bottomtransect.x = bathy.xgrid(indsecplt, :);
%
bottomtransect.depth = bathy.depth(indsecplt, :);


%% Calculate bottom slope

%
bottomtransect.dhdx = ...
                (bottomtransect.depth(3:end) - bottomtransect.depth(1:end-2)) ./ ...
                (bottomtransect.x(3:end) - bottomtransect.x(1:end-2));

%
bottomtransect.dhdx = bottomtransect.dhdx ./ 1000;

%
bottomtransect.dhdx = [NaN, bottomtransect.dhdx, NaN];


%% Interpolate stratification to the bottom

%
bottomtransect.N2 = interp1(profN2.z, profN2.N2, bottomtransect.depth);


%% Calculate internal wave characteristic

%
wvfreq = 2*pi/(3600*12.42);

%
lat0 = 41;
f0 = gsw_f(lat0);

%
bottomtransect.IWchar = sqrt((wvfreq^2 - f0^2) ./ ...
                             (bottomtransect.N2 - wvfreq^2));

%
bottomtransect.slopecritic = abs(bottomtransect.dhdx) ./ bottomtransect.IWchar;


%%
% --------------------------------------------------------------
% -------------------- NOW DO WKB-STRETCHING -------------------
% --------------------------------------------------------------


%%
%
bottomdepth_max = 4270;

%%

%
lN2inlims = (profN2.z >= 0) & (profN2.z <= bottomdepth_max);

%
N20ref = mean(profN2.N2(lN2inlims));


%% WKB-stretched bathymetry

%
bottomtransect.wkbstretched.N20 = N20ref;

%
bottomtransect.wkbstretched.depth = NaN(1, length(bottomtransect.depth));

%
for i = 1:length(bottomtransect.depth)
    
	%
    linabovebot_aux = (profN2.z <= bottomtransect.depth(i));
    
    %
    bottomtransect.wkbstretched.depth(i) = sum(sqrt(profN2.N2(linabovebot_aux))./sqrt(N20ref)) * profN2.dz;
    
end


%% Criticality of stretched bathymetry

%
bottomtransect.wkbstretched.dhdx = ...
                (bottomtransect.wkbstretched.depth(3:end) - bottomtransect.wkbstretched.depth(1:end-2)) ./ ...
                (bottomtransect.x(3:end) - bottomtransect.x(1:end-2));

%
bottomtransect.wkbstretched.dhdx = bottomtransect.wkbstretched.dhdx ./ 1000;

%
bottomtransect.wkbstretched.dhdx = [NaN, bottomtransect.wkbstretched.dhdx, NaN];

%
bottomtransect.wkbstretched.IWchar = sqrt((wvfreq^2 - f0^2) ./ ...
                                          (bottomtransect.wkbstretched.N20 - wvfreq^2));

%
bottomtransect.wkbstretched.slopecritic = ...
                                        abs(bottomtransect.wkbstretched.dhdx) ./ ...
                                        bottomtransect.IWchar;
                            
                                    
%%
% --------------------------------------------------------------
% ----- NOW RUN SCATTERING FOR AN INCIDENT MODE-1 INTERNAL -----
% --------------------------------------------------------------

%%

%
ridgeandflat = [3300, 850];


%% Basic definitions

% ------------------------------------------
%
H = ridgeandflat(1);

%
h0 = ridgeandflat(2) + 1;

%
gamma_cte = (H - h0)/H;

% Notice that after equation (A16), Klymak et al. (2013) makes
% a remark about gamma and singularities (I'm not sure what
% integer division of 1 means).

% ------------------------------------------
%
lat0 = T1loc.latitude;
f0 = gsw_f(lat0);

%
N = sqrt(bottomtransect.wkbstretched.N20);


% ------------------------------------------
% Internal modes
% Nmax = 80;
Nmax = 200;

%
nvec = 1:Nmax;
nvec = nvec(:);

%
mvec = 0:(Nmax-1);
%
mvec = mvec(:);


%% Forcing

% Magnitude of the incident internal tides
dn = zeros(Nmax, 1);

% % %
% % dn(1) = 2.705e-2;    % this magnitude is equal to the peak velocity in m/s
% %                      % (because the normal mode is defined by cos only)
% % dn(1) = 0.0191186206;    % energy flux of 0.5 kW/m
% dn(1) = 2.7038e-2;    % energy flux of 1.0 kW/m
% % dn(1) = 0.03022919;    % energy flux of 1.25 kW/m
% % dn(1) = 0.0331144;    % energy flux of 1.5 kW/m
% % dn(1) = 3.8237e-2;    % energy flux of 2 kW/m
dn(1) = 4.27505e-2;    % energy flux of 2.5 kW/m
% % dn(1) = 0.0448371396;    % energy flux of 2.75 kW/m
% % dn(1) = 0.0468308651;    % energy flux of 3 kW/m
% % dn(1) = 0.0505831155;    % energy flux of 3.5 kW/m
% % dn(1) = 0.0540756251;    % energy flux of 4.0 kW/m


%
U0 = 0;    % 0 surface tide forcing

% % % Surface tide forcing only (comment out the above)
% % U0 = 0.03;


% % Mmode = 38;    % 0.5
% % Mmode = 31;    % 1.0
% % Mmode = 29;    % 1.25
% % Mmode = 25;    % 2.0
Mmode = 23;    % 2.5
% % Mmode = 22;    % 3.0
% % Mmode = 20;    % 4.0



% % % % % Fluxes of [0.5 / 1 / 1.25 / 2 / 2.5 / 3 / 4]
Uscale = [0.0191186206, 2.7038e-2, 0.03022919, 0.0382, 0.0428, 0.0468308651, 0.0540756251];    % velocity amplitude of the mode-1
Dscale = [5.1143, 12.3670, 17.0624, 31.7843, 42.4887, 55.59203, 76.53222];
%
% Um to determine cut-off
Ucutoffscale = [0.059150857014622, ...    % 0.5
                0.075291189604048, ...    % 1.0
                0.079735300490465, ...    % 1.25
                0.093638033336456, ...    % 2.0
                0.100252598798004, ...    % 2.5
                0.109231746389776, ...    % 3.0
                0.117332362387293];       % 4.0


% % total_Tflux_arrested + total_Rflux_arrested
% % max(u_T_ridge_avg)

%% Check energy flux of the forcing

%
g_factor = 1030 * sqrt((N^2 - wvfreq^2)*(wvfreq^2 - f0^2)) / wvfreq;

% Depth-integrated flux (in W/m)
IncITflux = (H^2 ./ (pi*nvec)) .* (g_factor) .* (abs(dn).^2)./4;

IncITflux(1)

%% c vector

%
cm = (U0./mvec) .* sin((pi*gamma_cte) * mvec);

%
cm(1) = -U0 * pi * (1 - gamma_cte);


%%

%
[ngrid, mgrid] = meshgrid(nvec, mvec);


%% Get indices of the singular elements in the matrices.
% Singular locations are along the -1 diagonal
% (the one just below the main diagonal)

%
ble = eye(Nmax);
ble = [zeros(1, Nmax); ble(1:end-1, :)];
%
indsingular = find(ble);

%
nvec_sing = nvec(1:end-1);


%%

%
Amn = (ngrid .* sin((pi * gamma_cte) .* ngrid) .* cos((pi * gamma_cte) .* mgrid)) - ...
      (mgrid .* cos((pi * gamma_cte) .* ngrid) .* sin((pi * gamma_cte) .* mgrid));

%
Amn = Amn ./ (mgrid.^2 - ngrid.^2);

% Singular terms
Amn(indsingular) = ((pi * (1 - gamma_cte)) .* nvec_sing) - ...
                   (sin((pi * gamma_cte) .* nvec_sing) .* cos((pi * gamma_cte) .* nvec_sing));
Amn(indsingular) = Amn(indsingular) ./ (2 .* nvec_sing);


%%

%
Bmn = ngrid - ...
      (ngrid .* cos((pi * gamma_cte) .* ngrid) .* cos((pi * gamma_cte) .* mgrid)) - ...
      (mgrid .* sin((pi * gamma_cte) .* ngrid) .* sin((pi * gamma_cte) .* mgrid));

%
Bmn = Bmn ./ (mgrid.^2 - ngrid.^2);

% Singular terms
Bmn(indsingular) = - (sin((pi * gamma_cte) .* nvec_sing).^2);
Bmn(indsingular) = Bmn(indsingular) ./ (2.*nvec_sing);


%%  Solve for the internal-wave amplitudes

% ----------------------------------------------------
% bn are the amplitudes of the reflected waves
%
vector_RHS = cm - (Amn*dn);

%
full_matrix = Amn + Bmn;

%
bn = full_matrix \ vector_RHS;

% ----------------------------------------------------
% Now calculate amplitudes of the transmitted waves
% (see paragraph below BC's (A8)-(A10)

%
an = bn + dn;


%% Calculate flux of the transmitted waves

% These expressions should have the same form:
% % IncITflux = (H^2 ./ (pi*nvec)) .* (g_factor) .* (abs(dn).^2)./4;

% Depth-integrated flux (in W/m)
TransmitITflux = (H^2 ./ (pi*nvec)) .* (g_factor) .* (abs(an).^2)./4;


% Reflected waves
ReflectedITflux = (H^2 ./ (pi*nvec)) .* (g_factor) .* (abs(bn).^2)./4;



%%
% --------------------------------------------------------------
% --------------- NOW RECONSTRUCT VELOCITY FIELDS --------------
% --------------------------------------------------------------

%%

%
Xmax = 50 * 1000;
dx = 100;

%
xvec = 0 : dx : Xmax;
xvec = [-fliplr(xvec(2:end)), xvec];


%%

%
dz = H/200;

%
zvec = 0 : dz : H;
zvec = zvec(:);

%%

%
[xgrid, zgrid] = meshgrid(xvec, zvec);


%%

%
IWchar = sqrt( (wvfreq^2 - f0^2) ./ (N^2 - wvfreq^2) );

%
VertWaveNumber_vec = nvec .* (pi/H);
%
HorzWaveNumber_vec = IWchar * VertWaveNumber_vec;


%% Calculate eigenspeeds

%
cp_all = wvfreq ./ HorzWaveNumber_vec;

%
ce_all = cp_all .* sqrt(wvfreq^2 - f0^2)./ wvfreq;


%%

%
ntime = 1;
%
tvec = 0;

%%

xneg = xvec(xvec < 0);
xpos = xvec(xvec > 0);

%
[xneg_grid, ~] = meshgrid(xneg, zvec);
[xpos_grid, zgrid] = meshgrid(xpos, zvec);


%%

%
u_I = NaN([length(zvec), length(xpos), length(tvec), Nmax]);
u_R = NaN([length(zvec), length(xpos), length(tvec), Nmax]);
u_T = NaN([length(zvec), length(xneg), length(tvec), Nmax]);

%
tic
for i1 = 1:Nmax
    
    %
    vert_struct_aux = cos((pi*nvec(i1)/H) .* zgrid);
    
	%
    for i2 = 1:length(tvec)
        
        %
        wavy_osc_I_aux =  exp(1i .* ((HorzWaveNumber_vec(i1).*xpos_grid) + wvfreq*tvec(i2)));
        wavy_osc_R_aux =  exp(1i .* ((-HorzWaveNumber_vec(i1).*xpos_grid) + wvfreq*tvec(i2)));
        wavy_osc_T_aux =  exp(1i .* ((HorzWaveNumber_vec(i1).*xneg_grid) + wvfreq*tvec(i2)));
        
        %
        u_I(:, :, i2, i1) = dn(i1) .* vert_struct_aux .* wavy_osc_I_aux;
        %
        u_R(:, :, i2, i1) = bn(i1) .* vert_struct_aux .* wavy_osc_R_aux;
        %
        u_T(:, :, i2, i1) = an(i1) .* vert_struct_aux .* wavy_osc_T_aux;
        
        
    end    
end
toc

%
u_I = real(u_I);
u_R = real(u_R);
u_T = real(u_T);

%
disp('------------- DONE WITH VELOCITY RECONSTRUCTION -------------')


%%
% --------------------------------------------------------------
% ----------- NOW CALCULATE FIND THE ARRESTED MODES ------------
% ------------ (SECTION 4 IN THE PAPER, AND FIG. 4) ------------
% --------------------------------------------------------------

%% First recalculate the transmitted wave response, just where
% I want it (near the ridge) so I can calculate for "all" times.

%
xneg_aux = linspace(-1000, -10, 5);
%
[xneg_grid, zgrid] = meshgrid(xneg_aux, zvec);

%
tvec_aux = linspace(0, (2*pi/wvfreq), 20);

%
u_T_ridge = NaN([length(zvec), length(xneg_aux), length(tvec_aux), Nmax]);

%
tic
for i1 = 1:Nmax
    
    %
    vert_struct_aux = cos((pi*nvec(i1)/H) .* zgrid);
    
	%
    for i2 = 1:length(tvec_aux)
        
        %
        wavy_osc_T_aux =  exp(1i .* ((HorzWaveNumber_vec(i1).*xneg_grid) + wvfreq*tvec_aux(i2)));

        %
        u_T_ridge(:, :, i2, i1) = an(i1) .* vert_struct_aux .* wavy_osc_T_aux;
        
    end    
end
toc

%
u_T_ridge = real(u_T_ridge);
%
disp('-------- DONE WITH VELOCITY RECONSTRUCTION NEAR THE RIDGE --------')


%%

%
u_T_ridge_partial = sum(u_T_ridge(:, :, :, 1:Mmode), 4);


%% Average the velocity above the ridge

% Select the one closest to the top of the ridge
indcol = length(xneg_aux);

% Half the vertical wavelength of the critical mode
% (the vertical wavelength is 2*H/n / half is H/n)
dist_vert_aux = H / (Mmode);

%
ridge_depth = ridgeandflat(1) - ridgeandflat(2);

%
linlims_aux = (zvec < ridge_depth) & (zvec > (ridge_depth - dist_vert_aux));

%
u_T_ridge_prof = u_T_ridge_partial(:, indcol, :);
u_T_ridge_prof = squeeze(u_T_ridge_prof);

% Average over a half vertical wavelength
u_T_ridge_avg = mean(u_T_ridge_prof(linlims_aux, :), 1);

    


%%
% --------------------------------------------------------------
% --------------------------- FIGURES --------------------------
% --------------------------------------------------------------

%% 

% Integrate flux on each side of the mode cutoff
lhigharrestedmodes_aux = (nvec >= Mmode);

%
total_Tflux_arrested = sum(TransmitITflux(lhigharrestedmodes_aux));
total_Tflux_radiated = sum(TransmitITflux(~lhigharrestedmodes_aux));
%
total_Rflux_arrested = sum(ReflectedITflux(lhigharrestedmodes_aux));
total_Rflux_radiated = sum(ReflectedITflux(~lhigharrestedmodes_aux));


%
c_eigen_vec = sqrt(N20ref)*H./(nvec.*pi);
c_phasespeed_vec = c_eigen_vec .* (wvfreq) ./ sqrt(wvfreq^2 - f0^2);


%%

%
lnWd = 2;

%
cmaprb = load(fullfile(paper_directory(), 'figures', 'colormaps', ...
                                          'colormap_velocity.mat'));
cmaprb = cmaprb.cmapredblue;


%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.25, 0.25, 0.175, 0.6];
%
haxs = makeSubPlots(0.25, 0.2, 0.1, ...
                    0.02, 0.08, 0.075, 1, 3);
hold(haxs, 'on')

    % --------------------------------------------------
    %
    plot(haxs(1), bottomtransect.x, bottomtransect.depth, 'LineWidth', lnWd)
    plot(haxs(1), bottomtransect.x, bottomtransect.wkbstretched.depth, 'LineWidth', lnWd)

    %
    scatSZ = 25;
    %
    scatter(haxs(1), bottomtransect.x, bottomtransect.depth, scatSZ, bottomtransect.slopecritic, 'filled')
    scatter(haxs(1), bottomtransect.x, bottomtransect.wkbstretched.depth, scatSZ, bottomtransect.wkbstretched.slopecritic, 'filled')

	% Plot schematic
    overlayline(haxs(1), 'h', ridgeandflat(1), '-k', 'LineWidth', 4)
	%
    plot(haxs(1), [0, 0], [ridgeandflat(1), (ridgeandflat(1) - ridgeandflat(2))], 'k', 'LineWidth', 4)
    
    
    %
    colormap(haxs(1), cmaprb)
    
    %
    caxis(haxs(1), [0, 2])
        
	%
    hcb = colorbar(haxs(1));
        hcb.Position = [0.81, haxs(1).Position(2), 0.02, haxs(1).Position(4)];
        hcb.Label.String = '$\gamma$';
        hcb.Label.Interpreter = 'Latex';
        hcb.Label.FontSize = 40;
        hcb.Label.Rotation = 0;
        hcb.Label.VerticalAlignment = 'middle';
        hcb.FontSize = 14;
        hcb.Ticks = 0:0.5:2;
        hcb.Label.FontSize = 22;
        


    % --------------------------------------------------

    %
    u_T_ridge_allmodes = sum(u_T_ridge(:, :, :, :), 4);

    %
    [~, ind_max_aux] = max(u_T_ridge_avg);

    %
    plot(haxs(2), squeeze(u_T_ridge_allmodes(:, end, ind_max_aux)), zvec, ...
         'Color', 0.5.*[1, 1, 1])

    %
    plot(haxs(2), squeeze(u_T_ridge_partial(:, end, ind_max_aux)), zvec, '-k', 'LineWidth', 2)
    plot(haxs(2), [0, 0], [H, H-h0], 'k', 'LineWidth', 4)

    %
    plot(haxs(2), max(u_T_ridge_avg).*[1, 1], [ridge_depth, (ridge_depth-dist_vert_aux)], '-m', 'LineWidth', 2)
    
    %
    set(haxs(2), 'XLim', [-5, 5])

           
	%
    hline_aux = overlayline(haxs(2), 'v', c_eigen_vec(Mmode), '-', 'LineWidth', 2);
    ht_aux = text(haxs(2), c_eigen_vec(Mmode) - 0.03125, 1800, ['c$_{' num2str(Mmode) '}$']);
        ht_aux.FontSize = 22;
        ht_aux.Interpreter = 'Latex';
        ht_aux.Color = hline_aux.Color;
    %
    ht_aux_2 = text(haxs(2), 0.085, 2750, ['$U_{' num2str(Mmode) '}$']);
        ht_aux_2.FontSize = 22;
        ht_aux_2.Interpreter = 'Latex';
        ht_aux_2.Color = 'm';
	%
    overlayline(haxs(2), 'v', c_phasespeed_vec(Mmode), '--', 'LineWidth', 2)
    overlayline(haxs(2), 'h', ridgeandflat(1), '-k', 'LineWidth', 4)


    % --------------------------------------------------
    %
    plot(haxs(3), nvec, TransmitITflux + ReflectedITflux, '.-k')

    %
    ht_M_aux = text(haxs(3), 24, 250, ['$M = ' num2str(Mmode) '$']);
        ht_M_aux.FontSize = 18;
        ht_M_aux.Interpreter = 'Latex';
        ht_M_aux.Color = hline_aux.Color;
    
    %
    set(haxs(3), 'YLim', [1e-10, 1e10])
    overlayline(haxs(3), 'v', Mmode, '-', 'LineWidth', 2, 'Color', ht_aux.Color)
    

    % --------------------------------------------------
    %
    set(haxs, 'FontSize', 12, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs(1), 'XGrid', 'on', 'YGrid', 'on', 'YDir', 'reverse', 'XLim', [-31, 52], 'YLim', [0, 4300])
    %
    set(haxs(2), 'XLim', [-0.075, 0.125], 'YDir', 'reverse', 'YLim', [1400, (H+200)])
    %
    set(haxs(3), 'XLim', [7e-1, nvec(end)+10], 'XScale', 'log', 'XTick', [1, 10, 100], ...
                 'YLim', [1e-4, 3000], 'YScale', 'log', 'YTick', [1e-3, 1e-1, 10, 1e3])
    
    % --------------------------------------------------
    %
    lblFS = 14;
    %
    hxlb1 = xlabel(haxs(1), '[km]', 'Interpreter', 'Latex', 'FontSize', lblFS);
    ylabel(haxs(1), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
    %
    hxlb2 = xlabel(haxs(2), '[m s$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', lblFS);
    ylabel(haxs(2), 'Depth [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
    %
    xlabel(haxs(3), 'Mode number', 'Interpreter', 'Latex', 'FontSize', lblFS)
    ylabel(haxs(3), '[W m$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', lblFS)

    %
    hxlb1.Position(2) = 4850;
    hxlb2.Position(2) = 3790;
    
    %
    set(haxs, 'Color', 0.9.*[1, 1, 1])


    % --------------------------------------------------
    % Add letter labels
    hsqr_1 = fill(haxs(1), [-31, -31, -16, -16, -31], [4300, 3450, 3450, 4300, 4300], 'w');
    hsqr_2 = fill(haxs(2), [-0.075, -0.075, -0.038, -0.038, -0.075], [3500, 3100, 3100, 3500, 3500], 'w');
    hsqr_3 = fill(haxs(3), [0.7, 0.7, 2, 2, 0.7], [1e-4, 2.5e-3, 2.5e-3, 1e-4, 1e-4], 'w');
    %
    txtFS = 22;
    text(haxs(1), -24, 3920, 'a)', 'FontSize', txtFS, 'Interpreter', 'Latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    text(haxs(2), -0.0575, 3325, 'b)', 'FontSize', txtFS, 'Interpreter', 'Latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    text(haxs(3), 1.15, 4e-4, 'c)', 'FontSize', txtFS, 'Interpreter', 'Latex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')

    % ------------------
    htxt_1 = text(haxs(1), 18, 400, 'continental slope', 'FontSize', 10.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    htxt_2 = text(haxs(1), 34, 1500, {'bathymetry in';'scattering model'}, 'FontSize', 10.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    htxt_3 = text(haxs(1), 5, 3900, {'WKB-stretched';'bathymetry'}, 'FontSize', 10.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');


% -------------------------------------------------------
% -------------------------------------------------------
% Plot annotation arrows

%
haxs_arrows = axes('Position', haxs(1).Position);
hold(haxs_arrows, 'on')

%
shaftwidth = 0.14;
arrwhdwidth = 0.13;
arrwhdlen = 0.35;
%
xcoordbase = [0, 0, (1-arrwhdlen), (1-arrwhdlen), 1, (1-arrwhdlen), (1-arrwhdlen), 0];
ycoordbase = [-shaftwidth/2, +shaftwidth/2, +shaftwidth/2, +shaftwidth/2 + arrwhdwidth, 0, (-shaftwidth/2 - arrwhdwidth), -shaftwidth/2, -shaftwidth/2];

% --------------------------------
%
% axslength = 0.35;
% % axs_arrw_1 = axes('Position', [0.05, 0.7, axslength, axslength]);

    % Arrow 1
    %
    xy_coordlims = 0.475.*[xcoordbase(:).'; ycoordbase(:).'];
    rot_ang = 90;
    rot_matrix = [cosd(rot_ang), -sind(rot_ang); sind(rot_ang), cosd(rot_ang)];
    %
    xy_rotcoord = rot_matrix * xy_coordlims;
    %
    xy_rotcoord(1, :) = xy_rotcoord(1, :) - 0.6;
    xy_rotcoord(2, :) = xy_rotcoord(2, :) - 0.625;
    
    %
    hfill_ar1 = fill(haxs_arrows, xy_rotcoord(1, :), xy_rotcoord(2, :), 0.5.*[1, 1, 1]);
    hfill_ar1.LineStyle = 'none';

    
    % Arrow 2
    %
    xy_coordlims = [(0.45.*xcoordbase(:).'); (0.45*ycoordbase(:).')];
    rot_ang = -125;
    rot_matrix = [cosd(rot_ang), -sind(rot_ang); sind(rot_ang), cosd(rot_ang)];
    %
    xy_rotcoord = rot_matrix * xy_coordlims;
    %
    xy_rotcoord(1, :) = xy_rotcoord(1, :) - 0.3;
    xy_rotcoord(2, :) = xy_rotcoord(2, :) + 0.7;
    
    %
    hfill_ar2 = fill(haxs_arrows, xy_rotcoord(1, :), xy_rotcoord(2, :), 0.5.*[1, 1, 1]);
    hfill_ar2.LineStyle = 'none';
    
    
    % Arrow 3
    %
    xy_coordlims = [(0.6.*xcoordbase(:).'); (0.475*ycoordbase(:).')];
    rot_ang = -90;
    rot_matrix = [cosd(rot_ang), -sind(rot_ang); sind(rot_ang), cosd(rot_ang)];
    %
    xy_rotcoord = rot_matrix * xy_coordlims;
    %
    xy_rotcoord(1, :) = xy_rotcoord(1, :) + 0.6;
    xy_rotcoord(2, :) = xy_rotcoord(2, :) + 0.1;
    
    %
    hfill_ar3 = fill(haxs_arrows, xy_rotcoord(1, :), xy_rotcoord(2, :), [0, 0, 0]);
    hfill_ar3.LineStyle = 'none';

%
set(haxs_arrows, 'DataAspectRatio', [1, 1, 1])
set(haxs_arrows, 'XLim', [-1, 1], 'YLim', [-1, 1])

%
haxs_arrows.Visible = 'off';
    
%
exportgraphics(gcf, 'figure09.pdf', 'Resolution', 300)
