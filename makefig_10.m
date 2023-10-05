%% Make velocity vs. dissipation scaling figure


clear
close all


%% Load results

load(fullfile(paper_directory(), 'data', 'mooring', 'level_3', 'T1_dissipation_percycle.mat'))


%% Select data

%
time_break_refs = [26.5, 55];

%
lintime_1 = (phasecycles.scales.timemidref > 18.5) & ...
            (phasecycles.scales.timemidref < time_break_refs(1));
%
lintime_2 = (phasecycles.scales.timemidref > time_break_refs(1)) & ...
            (phasecycles.scales.timemidref < time_break_refs(2));
%
lintime_3 = (phasecycles.scales.timemidref > time_break_refs(2)) & ...
             phasecycles.scales.timemidref < 60.5;


%
x_fit_data_2a = phasecycles.scales.U(lintime_2);
x_fit_data_2b = phasecycles.scales.V(lintime_2);

%
x_fit_data_1a = phasecycles.scales.U(lintime_1);
x_fit_data_1b = phasecycles.scales.V(lintime_1);

%
x_fit_data_3a = phasecycles.scales.U(lintime_3);
x_fit_data_3b = phasecycles.scales.V(lintime_3);

%
y_fit_data_1 = phasecycles.scales.dissipation(lintime_1);
y_fit_data_2 = phasecycles.scales.dissipation(lintime_2);
y_fit_data_3 = phasecycles.scales.dissipation(lintime_3);


%% Compute statistics using only U

% -----------------------------------------
% First period
x_fit_var_1 = log10(x_fit_data_1a(:));
y_fit_var_1 = log10(y_fit_data_1(:));
%
x_fit_matrix = [ones(length(x_fit_var_1), 1), x_fit_var_1];
%
[line_params_fit_1a, params_err_1a] = regress(y_fit_var_1, x_fit_matrix);

% -----------------------------------------
% Second period
x_fit_var_2 = log10(x_fit_data_2a(:));
y_fit_var_2 = log10(y_fit_data_2(:));
%
x_fit_matrix_2 = [ones(length(x_fit_var_2), 1), x_fit_var_2];
%
[line_params_fit_2a, params_err_2a] = regress(y_fit_var_2, x_fit_matrix_2);

%
rfit_2a = corr(x_fit_var_2(:), y_fit_var_2(:));


% -----------------------------------------
% Third period
x_fit_var_3 = log10(x_fit_data_3a(:));
y_fit_var_3 = log10(y_fit_data_3(:));
%
x_fit_matrix_3 = [ones(length(x_fit_var_3), 1), x_fit_var_3];
%
[line_params_fit_3a, params_err_3a] = regress(y_fit_var_3, x_fit_matrix_3);

%% Compute statistics using UV

% -----------------------------------------
% First period
x_fit_var_1 = log10(x_fit_data_1b(:));
y_fit_var_1 = log10(y_fit_data_1(:));
%
x_fit_matrix = [ones(length(x_fit_var_1), 1), x_fit_var_1];
%
[line_params_fit_1b, params_err_1b] = regress(y_fit_var_1, x_fit_matrix);

% -----------------------------------------
% Second period
x_fit_var_2 = log10(x_fit_data_2b(:));
y_fit_var_2 = log10(y_fit_data_2(:));
%
x_fit_matrix_2 = [ones(length(x_fit_var_2), 1), x_fit_var_2];
%
[line_params_fit_2b, params_err_2b] = regress(y_fit_var_2, x_fit_matrix_2);

%
rfit_2b = corr(x_fit_var_2(:), y_fit_var_2(:));

% -----------------------------------------
% Third period
x_fit_var_3 = log10(x_fit_data_3b(:));
y_fit_var_3 = log10(y_fit_data_3(:));
%
x_fit_matrix_3 = [ones(length(x_fit_var_3), 1), x_fit_var_3];
%
[line_params_fit_3b, params_err_3b] = regress(y_fit_var_3, x_fit_matrix_3);


%% Compute statistics for the full timeseries


% --------------------------------------------
%
x_fit_var_all_a = log10([x_fit_data_1a(:); x_fit_data_2a(:); x_fit_data_3a(:)]);
y_fit_var_all_a = log10([y_fit_data_1(:); y_fit_data_2(:); y_fit_data_3(:)]);
%
x_fit_matrix_all_a = [ones(length(x_fit_var_all_a), 1), x_fit_var_all_a];
%
[line_params_fit_all_a, params_err_all_a] = regress(y_fit_var_all_a, x_fit_matrix_all_a);


%
rfit_all_a = corr(x_fit_var_all_a(:), y_fit_var_all_a(:));


% --------------------------------------------
%
x_fit_var_all_b = log10([x_fit_data_1b(:); x_fit_data_2b(:); x_fit_data_3b(:)]);
y_fit_var_all_b = log10([y_fit_data_1(:); y_fit_data_2(:); y_fit_data_3(:)]);
%
x_fit_matrix_all_b = [ones(length(x_fit_var_all_b), 1), x_fit_var_all_b];
%
[line_params_fit_all_b, params_err_all_b] = regress(y_fit_var_all_b, x_fit_matrix_all_b);


%
rfit_all_b = corr(x_fit_var_all_b(:), y_fit_var_all_b(:));



%% Make the figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.4105    0.1258    0.2141    0.5583];
%
haxs = makeSubPlots(0.25, 0.25, 0.025, ...
                    0.05, 0.15, 0.15, 1, 2);
hold(haxs, 'on')


    % ------------------------------------------------
    % Plot each period

    %
    mkSZdots = 20;

    % ------------------------------
    %
    plot(haxs(1), x_fit_data_1a, y_fit_data_1, '.r', 'MarkerSize', mkSZdots);
    plot(haxs(1), x_fit_data_3a, y_fit_data_3, '.b', 'MarkerSize', mkSZdots);
    %
    plot(haxs(1), x_fit_data_2a, y_fit_data_2, '.k', 'MarkerSize', mkSZdots);

    %
    Upick = linspace(1e-2, 1.1e-1, 50);
    hp_KnifeEdge = plot(haxs(1), Upick, (1/5000)*(10^5.2339)*Upick.^(2.634728), '--g', 'LineWidth', 2);

    % ------------------------------
    %
    hplt_1 = plot(haxs(2), x_fit_data_1b, y_fit_data_1, '.r', 'MarkerSize', mkSZdots);
    hplt_2 = plot(haxs(2), x_fit_data_3b, y_fit_data_3, '.b', 'MarkerSize', mkSZdots);
    %
    hplt_3 = plot(haxs(2), x_fit_data_2b, y_fit_data_2, '.k', 'MarkerSize', mkSZdots);


    % ------------------------------------------------
    % Plot best fit line

    %
    xplt_aux = [0.025, 0.15];
    yplt_aux = 10.^(line_params_fit_2a(1)).*xplt_aux.^(line_params_fit_2a(2));

    %
    hp_bestfit = plot(haxs(1), xplt_aux, yplt_aux, '--k', 'LineWidth', 2);


    %
    xplt_aux = [0.025, 0.2];
    yplt_aux = 10.^(line_params_fit_2b(1)).*xplt_aux.^(line_params_fit_2b(2));


    %
    plot(haxs(2), xplt_aux, yplt_aux, '--k', 'LineWidth', 2)
        
% ------------------------------------------------
% Legend of the regression slope

%
txtFS = 14;

%
htxt_a_aux = text(haxs(1), 0.025, 5.8e-4, ['slope = ' num2str(line_params_fit_2a(2), '%.1f') '\pm' num2str(line_params_fit_2a(2) - params_err_2a(2, 1), '%.1f')]);                                
% % % % htxt_b_aux = text(haxs(1), 0.06, 0.775e-3, ['r^2 = ' num2str(rfit_2a.^2, '%.1f')]);
% % htxt_b_aux = text(haxs(1), 0.0625, 1.2e-3, ['r^2 = ' num2str(rfit_2a.^2, '%.1f')]);
    %
    htxt_a_aux.FontSize = txtFS;
    htxt_a_aux.Color = [0, 0, 0];
    %
% %     htxt_b_aux.FontSize = txtFS;
% %         htxt_b_aux.Interpreter = 'Latex';



%
htxt_c_aux = text(haxs(2), 1.02e-2, 5.8e-4, ['slope = ' num2str(line_params_fit_2b(2), '%.1f') ...
                                    '\pm' num2str(line_params_fit_2b(2) - params_err_2b(2, 1), '%.1f')]);                                
% % htxt_d_aux = text(haxs(2), 1.1e-2, 0.9e-3, ['r^2 = ' num2str(rfit_2b.^2, '%.1f')]);
    %
    htxt_c_aux.FontSize = txtFS;
% %     htxt_c_aux.Color = [0, 0, 0];
    %
% %     htxt_d_aux.FontSize = txtFS;
% %         htxt_b_aux.Interpreter = 'Latex';


        

% ------------------------------
% Plot legend
%
hleg = legend(haxs(2), [hplt_1, hplt_3, hplt_2, hp_KnifeEdge], ...
                'Period I', 'Period II', 'Period III', 'Knife edge', 'Location', 'NorthWest');
% %             hleg.Position = [0.2515, 0.82, 0.2215, 0.13];
    hleg.Position = [0.2488    0.3333    0.2994    0.1410];
    hleg.FontSize = 13;
    
    
% ------------------------------
%
set(haxs, 'FontSize', 16, 'Box', 'on', ...
                    'XGrid', 'on', 'YGrid', 'on', ...
                    'XScale', 'log', 'YScale', 'log', ...
                    'YLim', [5e-4, 1e-2])
set(haxs, 'XLim', [0.01, 0.15])
%
set(haxs, 'XTick', [1e-2, 1e-1, 1])
set(haxs, 'Color', 0.65.*[1, 1, 1])


%
lblFS = 18;
%
hxlbl_1 = xlabel(haxs(1), '[m s$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', lblFS);
hxlbl_2 = xlabel(haxs(2), '[m s$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', lblFS);

%
hylbl_1 = ylabel(haxs(1), '[W m${-^2}$]', 'Interpreter', 'Latex', 'FontSize', lblFS);
hylbl_2 = ylabel(haxs(2), '[W m${-^2}$]', 'Interpreter', 'Latex', 'FontSize', lblFS);

%
hxlbl_1.Position(2) = 3.6e-4;
hxlbl_2.Position(2) = 3.6e-4;
%
hylbl_1.Position(1) = 6e-3;
hylbl_2.Position(1) = 6e-3;

%
title(haxs(1), '$U_{\mathrm{D}_2}$ vs. dissipation', 'Interpreter', 'Latex', 'FontSize', 16)
title(haxs(2), '$V_{\mathrm{D}_2}$  vs. dissipation', 'Interpreter', 'Latex', 'FontSize', 16)


% ------------------------------
% Add letter labelling
%
text(haxs(1), 1.025e-2, 1.2e-2, 'a)', 'FontSize', 16)
text(haxs(2), 1.025e-2, 1.2e-2, 'b)', 'FontSize', 16)


%% Save figure

exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure10.pdf'), 'Resolution', 300)





