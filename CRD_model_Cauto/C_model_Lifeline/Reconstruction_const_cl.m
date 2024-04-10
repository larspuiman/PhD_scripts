clearvars
format compact
 close all

%% lifeline CFD data
lifeline_data = load('lifeline_Ac90_8tm_part_8.mat');
lifeline_data = lifeline_data.lifeline_data;
conc_data = lifeline_data(:,12:12+11);
dcdt_data = lifeline_data(:,12+16 : 12+16+11);
t_data = lifeline_data(:,3);
t_data = t_data - t_data(1);
NP = lifeline_data(:,24);
compounds = {'CO', 'H_2','CO_2', 'AcT', 'EtOH', 'BDO', 'AcCoA', 'NADH', ...
    'NADPH', 'Ac^-_{IC}', 'For^-', 'Fd_{red}^{2-}', 'X'};
loc_data = lifeline_data(:,4:6);
locs = {'x','y','z'};

tmax = 1000;
ind_max = find(t_data < tmax);
conc_data = conc_data(ind_max, :);
t_data = t_data(ind_max, :);
loc_data = loc_data(ind_max, :);
NP = NP(ind_max, :);
dcdt_data = dcdt_data(ind_max, :);
mean_CO = round(mean(conc_data(:,1)), 2);
smoothingfactor = 0;
conc_data_s = smoothdata(conc_data, 'movmean', 'SmoothingFactor', smoothingfactor);

%% identify low concentration zones of CO:
low_conc_zones = [0];
previous_low = 0;
conc_CO_threshold = mean_CO/2;
conc_CO_threshold = 0.025;
for i = 1:length(conc_data(:,1))
    if conc_data_s(i,1) < conc_CO_threshold
        if (previous_low == 0)
            low_conc_zones = [low_conc_zones, t_data(i),t_data(i),];
        end
        previous_low = 1;
    elseif (conc_data_s(i,1) > conc_CO_threshold) && (previous_low == 1)
        low_conc_zones = [low_conc_zones, t_data(i),t_data(i)];
        previous_low = 0;
    else
        previous_low = 0;
    end
    
end
low_conc_zones(1) = [];
%% lifeline steady state
conc_data_ss = load('Const_cl_model/conc_data_cCO_0.29.dat');
dcdt_data_ss = load('Const_cl_model/dcdt_data_cCO_0.29.dat');
t_data_ss = load('Const_cl_model/t_data_cCO_0.29.dat');
conc_data_ss = conc_data_ss';
dcdt_data_ss = dcdt_data_ss';
t_data_ss = t_data_ss';


%% lifeline 1-way
conc_data_1way = load('Lifeline_reconstruct_C_1way/Results/conc_data_1way_ll8.dat');
dcdt_data_1way = load('Lifeline_reconstruct_C_1way/Results/dcdt_data_1way_ll8.dat');
t_data_1way = load('Lifeline_reconstruct_C_1way/Results/t_data_1way_ll8.dat');
conc_data_1way = conc_data_1way';
dcdt_data_1way = dcdt_data_1way';
t_data_1way = t_data_1way';






%% plot and compare
fig = figure(1);
set(fig,'DefaultAxesFontSize',16);
set(fig,'DefaultAxesFontName','times');
set(fig,'position', 1e3*[   -1.9190   -0.1750    1.9200    0.9648]);
for i = 1:size(conc_data,2)
    subplot(4,3,i)
    plot(t_data, conc_data_s(:,i), 'r-','LineWidth',1.5); hold on;
    plot(t_data_ss, conc_data_ss(:,i), 'b--'); hold on;
    plot(t_data_1way, conc_data_1way(:,i), 'k--'); hold on;
     xlim([0, 1000]);
     ylims = ylim();    
    if i == 8
        ylim([0, 0.06]);
    elseif i == 9
        ylim([0, 0.2]);
    end
    for ii = 0:(length(low_conc_zones)-1)/4
        fill([low_conc_zones(ii*4+1), low_conc_zones(ii*4+3), low_conc_zones(ii*4+4), low_conc_zones(ii*4+2)],...
            [ylims(1), ylims(1), ylims(2),ylims(2)],[0.5,0.5,0.5],'FaceAlpha',0.4,'EdgeColor','none');
    end
    ylabel(['{\itc}_{' compounds{i} '}'],'FontSize',15);
%         title(compounds{i});
    xlabel('{\itt}','FontSize',15)
%     set(gcf,'position', [117.8 122.6 1216 640.4]);
end
% saveas(gcf,'conc_ss_compare_reconstruction.png')
exportgraphics(gcf,'conc_ss_compare_reconstruction.png','Resolution',600);

