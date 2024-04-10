% this is the model of the scale-down simulator of the EL-GLR with 5 g/L
% it applies 0.05 vmm in the stirred tank. 

% The following results were obtained with conditions (table X in paper)
%   Peak:       rpm, kla, PV:   910 /min, 153 h-1, 2307 W/m3
%   Valley:     rpm, kla, PV:   20 /min, 1.129 h-1, 0.1086 W/m3
%   Startup:    rpm, kla, PV:   75 /min, 6.1733 h-1, 7.5901 W/m3
%   Recycle:    rate = 0.5, cx = 0.5414 g/L

clearvars
close all

% Do you want to plot and save figures and results? 1 = yes, 0 = no. 
plot_figs = 1;
save_figs = 0;
save_data = 0;

% 1) assign number of peaks to load - determines duration of simulation
Npeaks = 2000;

% 2) assign the parameters of the scale-down simulator in the function below
% "fun_CSTR" 

% 3) put in your favorite kinetic model in the CO_uptake and H2_uptake
% functions. The resulting q-rates for uptake should be positive-signed and in
% unit of mol/molx/h. 

% 4) put in your reactor dimensions and kLa formula in the "fKla" and
% "Volumes" functions 

% 5) Decide which industrial condition you want to replicate (5, 10 or 25
% gl) by putting adjusting the filenames in lines 128

% data path
DATA_PATH = "Underlying_data";
load(strcat(DATA_PATH,'/prob-peaks-val-clCO-5-gl'));
load(strcat(DATA_PATH,'/prob-peaks-val-clH2-5-gl'));

% run start-up ODE of reactor
t0 = 0; tend1 = 1500; dt = 0.1;        % h
tspan = t0:dt:tend1;                   % timespan (in hours)
MWx = 24.6;                            % molecular (bio)mass (g/mol)
c0 = [0.1, 0.03, 7.4, 0.05/MWx*1000];  % inlet concentrations [CO, H2, CO2, X] (mM)
y0 = [0.5, 0.2, 0.3];                  % gas mole fractions in dispersed gas phase [CO, H2, CO2]             
yh0 = [0.5, 0.2, 0.3];                 % gas mole fractions in headspace [CO, H2, CO2]      
x0 = [c0, y0, yh0];                    % initial values for ODE

%% start-up ODE system
[t, x_startup] = ode15s(@(t,x)fun_CSTR(t,x,[0,0]), tspan, x0);

% unpack results from start-up ODE
c = x_startup(:,1:4);                  % dissolved concentrations [CO, H2, CO2, X] (mM)       
y = x_startup(:,5:7);                  % gas mole fractions in dispersion [CO, H2, CO2]
yh = x_startup(:,8:10);                % gas mole fractions in headspace [CO, H2, CO2]      

if plot_figs == 1
    figure('DefaultAxesFontSize',12)
    % CO over time
    subplot(1,3,1)
    plot(t, c(:,1),'r-', 'LineWidth',2)
    xlabel('Time (h)')
    ylabel('c_{L,CO} (mM)')

    % H2 over time
    subplot(1,3,2)
    plot(t, c(:,2),'r-', 'LineWidth',2)
    xlabel('Time (h)')
    ylabel('c_{L,H_2} (mM)')
   
    % Biomass over time
    subplot(1,3,3)
    plot(t, c(:,4)*MWx/1000,'r-', 'LineWidth',2)
    xlabel('Time (h)')
    ylabel('c_{L,X} (g/L)')
    if save_figs == 1, saveas(gcf,'Startup_phase.png'); end
    if save_data == 1, save('Startup_data.mat','t', 'c'); end
end


%% load peaks and valleys

% create the variational peaks and valleys 
    % all times in seconds now!!
[t_new_peak, t_old_peak, valley_times, peak_times] = Peak_development(Npeaks, peak_time_bins_clCO, prob_peak_time_clCO, valley_time_bins_clCO, prob_valley_time_clCO); 

% specificy time vector based on the total duration
time_vec = 1:1:ceil(max(t_new_peak));         % time vector (in seconds!)
p = zeros(length(time_vec),1);                % Boolean to store whether we are in peak (1) or valley (0)

% fill in the booleans
for i = 1:length(time_vec)
    % if we are in a peak
    if any(time_vec(i) >= t_new_peak & time_vec(i) <= [t_old_peak, t_old_peak(end)])
        p(i) = 1; 
    % or we are in a valley
    else 
        p(i) = 0;
    end
end

% figure to illustrate peak/valley formation
if plot_figs == 1
    figure()
    subplot(1,3,1)
    % time in peaks 
    plot(1:Npeaks, peak_times, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','r')
    xlabel('Peak number'); ylabel('Peak duration (s)');
    ylim([0, 20]); grid on;

    subplot(1,3,2)
    % time valleys
    plot(1:(Npeaks+1), valley_times, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor','r')
    xlabel('Valley number'); ylabel('Time in valley (s)');
    ylim([0, 100]); grid on;

    subplot(1,3,3)
    plot(time_vec, p, 'k-')
    xlabel('Time (s)'); ylim([-0.2, 1.2]);
    yticks([0 1])
    yticklabels({'Valley','Peak'})
    set(gcf,'position', [ -1289.4  429.8   1176  293.6]);
    % saveas(gcf,'Peak_development.png')
end

%% check divergence from CFD data
KL_valley = KL_divergence(prob_valley_time_clCO, valley_times, valley_time_bins_clCO);
KL_peak = KL_divergence(prob_peak_time_clCO, peak_times, peak_time_bins_clCO);


%% Fluctuations with multiple ODEs
% get initial values from start-up
x0 = x_startup(end,:); 
t0 = 0; dt = 1/3600;                        % time starts at 0 s, steps of 1 s

% define timespan for each peak and valley
tvector_pk_val = zeros(1, (2*Npeaks+1));    % the total amount of time for all the the fluctuations
tvector_pk_val(1) =  valley_times(1);       % first we start in a valley
% fill in the timevector for the rest of the peaks and valleys
for i=2:(2*Npeaks+1)   
    if mod(i,2) == 1  % now uneven number so we add the valley time
        tvector_pk_val(i) = tvector_pk_val(i-1) + valley_times(i/2+0.5);
    elseif mod(i,2) == 0 % even number so we add the peak time
        tvector_pk_val(i) = tvector_pk_val(i-1) + peak_times(i/2);
    end
end

% allocate matrices for the output data
t_loop = zeros(1,round(tvector_pk_val(end))); 
c_loop = zeros(round(tvector_pk_val(end)),4);
y_loop = zeros(round(tvector_pk_val(end)),3);
yh_loop = zeros(round(tvector_pk_val(end)),3);
ind_start = 0;                  % counter for storing data

% run ODE for each peak or valley
for i = 1:(2*Npeaks+1)
    % define peak/valley time-span (in hours!, because thats the unit of the
    % ODE solver)
    tend = round(tvector_pk_val(i))/3600;           % end of peak time in hours
    tspan = t0:dt:tend;                             % time span for ODE solver (hours!)
    
    % run ODE for a peak or valley - the p is an additional input to
    % indicate whether we are in peak / valley, stirrer speed is adjusted
    % accordingly
    [t, x] = ode15s(@(t,x)fun_CSTR(t,x,p), tspan, x0);   
    
    % unpack results
    c = x(:,1:4);                           % dissolved concentrations [CO, H2, CO2, X] (mM)       
    y = x(:,5:7);                           % gas mole fractions in dispersion [CO, H2, CO2]
    yh = x(:,8:10);                         % gas mole fractions in headspace [CO, H2, CO2]      
    
    % store t c and y values
    ind_end = length(t)+ind_start;
    t_loop(ind_start+1:ind_end) = t;
    c_loop(ind_start+1:ind_end,:) = c;    
    y_loop(ind_start+1:ind_end,:) = y;    
    yh_loop(ind_start+1:ind_end,:) = yh;    

    % update time and concentrations for next peak/valley
    t0 = t(end);
    x0 = x(end,:);
    ind_start = ind_end;
end
    
% % retrieve uptake rates
% q = zeros(size(c_loop)); pode = zeros(size(t_loop));
% for i = 1:length(t_loop)
%     [dcdt, q(i,:), pode(i)] = fun(t_loop(i), c_loop(i,:), time_pks, p);
% end

%% Plot ODE results
if plot_figs == 1
    figure()
    % plot dissolved CO concentrations
    subplot(1,3,1)
    plot(t_loop*3600, c_loop(:,1),'r-', 'LineWidth',2)
    xlabel('Time (s)')
    ylabel('c_{L,CO} (mM)');
    xlim([0, max(t_loop)*3600])

    subplot(1,3,2)
    % plot CO gas fraction in dispersion
    plot(t_loop*3600, y_loop(:,1),'r-', 'LineWidth',2)
    xlabel('Time (s)')
    ylabel('y_{CO}'); 
    xlim([0, max(t_loop)*3600])
    set(gcf,'position',[410.6000  342.6000  922.4000  420.0000]);
    
    subplot(1,3,3)
    % plot CO gas fraction in headspace
    plot(t_loop*3600, yh_loop(:,1),'r-', 'LineWidth',2)
    xlabel('Time (s)')
    ylabel('y_{H,CO}'); 
    xlim([0, max(t_loop)*3600])
    set(gcf,'position',[410.6000  342.6000  922.4000  420.0000]);
end


%% Analyze the scale-down lifelines for CO and H2

factor = 1;                             % discrimination factor to distinguish between peak and valley
% co performance
conc_SD = c_loop(:,1);
dt = (t_loop(2)-t_loop(1))*3600;
[tCO_pks, tCO_val, cCO_pks, cCO_val] = lifeline_analysis(conc_SD, factor, Npeaks, dt);

% compare CO performance with CFD results using KL divergence
KL_cCO_peak = KL_divergence(prob_peak_conc_clCO, cCO_pks, peak_conc_bins_clCO);
KL_cCO_val = KL_divergence(prob_valley_conc_clCO, cCO_val, valley_conc_bins_clCO);
KL_tCO_peak = KL_divergence(prob_peak_time_clCO, tCO_pks, peak_time_bins_clCO);
KL_tCO_val = KL_divergence(prob_valley_time_clCO, tCO_val, valley_time_bins_clCO);

% H2 performance
conc_SD = c_loop(:,2);
[tH2_pks, tH2_val, cH2_pks, cH2_val] = lifeline_analysis(conc_SD, factor, Npeaks, dt);
KL_cH2_peak = KL_divergence(prob_peak_conc_clH2, cH2_pks, peak_conc_bins_clH2);
KL_cH2_val = KL_divergence(prob_valley_conc_clH2, cH2_val, valley_conc_bins_clH2);
KL_tH2_peak = KL_divergence(prob_peak_time_clH2, tH2_pks, peak_time_bins_clH2);
KL_tH2_val = KL_divergence(prob_valley_time_clH2, tH2_val, valley_time_bins_clH2);

%% plot  (histogram for SD data and line graph for CFD data
figure()
% plot for concentrations in peak (blue) and valley (red)
subplot(2,2,[1 2])
histogram(cCO_pks, peak_conc_bins_clCO, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.255, 0.412, 0.882]); hold on;
histogram(cCO_val, valley_conc_bins_clCO, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.545, 0,0]); 
plot(peak_conc_bins_clCO(1:end-1),prob_peak_conc_clCO, '-',...
    'LineWidth',2, 'Color',[0.255, 0.412, 0.882]); 
plot(valley_conc_bins_clCO(1:end-1),prob_valley_conc_clCO, '-', ...
    'LineWidth',2,'Color',[0.545, 0,0]); 
xlabel('CO concentration c_{L,CO}, (mol m^{-3})')
ylabel('probability')

% plot for time in valley
subplot(2,2,3)
histogram(tCO_val, valley_time_bins_clCO, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.545, 0,0]); hold on;
plot(valley_time_bins_clCO(1:end-1),prob_valley_time_clCO, '-', ...
        'LineWidth',2,'Color',[0.545, 0,0]); 
xlabel('Time in valley t_{valley,CO}, (s)')
ylabel('probability')

subplot(2,2,4)
% plot for time in peak
histogram(tCO_pks, peak_time_bins_clCO, 'Normalization','Probability', ...
    'LineStyle','none', 'FaceColor',[0.255, 0.412, 0.882]); hold on; hold on;
plot(peak_time_bins_clCO(1:end-1),prob_peak_time_clCO,'-',...
    'LineWidth',2, 'Color',[0.255, 0.412, 0.882]); 
xlabel('Time in peak t_{peak,CO}, (s)')
ylabel('probability')
if save_figs == 1; saveas(gcf,'CO_Scale-down-5gl.png'); end


figure() % same figure for H2
subplot(2,2,[1 2])
histogram(cH2_pks, peak_conc_bins_clH2, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.255, 0.412, 0.882]); hold on;
histogram(cH2_val, valley_conc_bins_clH2, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.545, 0,0]); 
plot(peak_conc_bins_clH2(1:end-1),prob_peak_conc_clH2, '-',...
    'LineWidth',2, 'Color',[0.255, 0.412, 0.882]); 
plot(valley_conc_bins_clH2(1:end-1),prob_valley_conc_clH2, '-', ...
    'LineWidth',2,'Color',[0.545, 0,0]); 
xlabel('H_2 concentration c_{L,H2}, (mol m^{-3})')
ylabel('probability')

subplot(2,2,3)
histogram(tH2_val, valley_time_bins_clH2, 'Normalization','Probability',...
    'LineStyle','none', 'FaceColor',[0.545, 0,0]); hold on;
plot(valley_time_bins_clH2(1:end-1),prob_valley_time_clH2, '-', ...
        'LineWidth',2,'Color',[0.545, 0,0]); 
xlabel('Time in valley t_{valley,H2}, (s)')
ylabel('probability')

subplot(2,2,4)
histogram(tH2_pks, peak_time_bins_clH2, 'Normalization','Probability', ...
    'LineStyle','none', 'FaceColor',[0.255, 0.412, 0.882]); hold on; hold on;
plot(peak_time_bins_clH2(1:end-1),prob_peak_time_clH2,'-',...
    'LineWidth',2, 'Color',[0.255, 0.412, 0.882]); 
xlabel('Time in peak t_{peak,H2}, (s)')
ylabel('probability')
if save_figs == 1;saveas(gcf,'H2_Scale-down-5gl'); end

%% save data
KL_all = [KL_cCO_peak, KL_cCO_val, KL_tCO_peak, KL_tCO_val, KL_cH2_peak, KL_cH2_val, KL_tH2_peak, KL_tH2_val];
if save_data == 1
    save('BenchScale-5-gl.mat','KL_all','cCO_pks','cCO_val', 'tCO_pks', 'tCO_val',...
         'cH2_pks','cH2_val', 'tH2_pks', 'tH2_val','c_loop','t_loop')
end
sum(KL_all)

%% Lifeline analysis
function [t_pks, t_val, c_pks, c_val] = lifeline_analysis(conc_SD, factor, Npeaks, dt)
% lifeline_analysis     to discriminate the peaks and valleys from a lifeline. 
%   *conc_SD* is a vector of which the peaks and valleys are to be
%   analyzed. Discrimination is done by imposing a *factor* that 
%   is based on the difference of the minimum/maximum 
%   concentration in peak/valley with the average concentration.  *Npeaks* is
%   the total amount of peaks imposed in the analysis, required for memory
%   allocation and *dt* is the time-step between each value in the
%   concentration vector. 
%   outputs the time and concentrations in the peaks and valleys

% retrieve average concentration during the peak/valley fluctuations
c_av = mean(conc_SD);                   % concentration (mM)

% check difference with average for concentrations
c_diff = conc_SD-c_av;                 % difference from average (mM)
%   The factor can be adjusted to sharpen discrimination between peaks and valleys
min_peak_c = c_av*factor;              % minimum concentration for a peak - mM
max_val_c = c_av/factor;               % maximum concentration for a valley - mM

% find indices for sign change 
%       (+1 = from negative to pos ==> valley to peak, 
%       -1 = from pos to negative ==> peak to valley.
signchange_ind = sign(diff(sign(c_diff)));

% allocate memory
c_sum = 0; t_sum = 0; counter = 0; c_max = min_peak_c; c_min = max_val_c;
c_val = zeros(1,Npeaks*5); t_val = zeros(1,Npeaks*5);
c_pks = zeros(1,Npeaks*5); t_pks = zeros(1,Npeaks*5);
ind_val = 1; ind_pks  = 1;
t_crit = 2;             % only store after 2 s in peak or valley

% analyze the lifeline by comparing the concentrations at a moment with the
% concentration min/max in a peak/valley
for i = 1:length(conc_SD)-1
    c_sum = c_sum + conc_SD(i);          % cumulative concentration so far in peak/valley
    t_sum = t_sum + dt;                  % cumulative concentration so far in peak/valley
    c_max = max([c_max, conc_SD(i)]);    % to check whether we are in a valley or going out of the valley
    c_min = min([c_min, conc_SD(i)]);    % to check whether we are in a peak or going out of the peak
    counter = counter + 1;               
    
    if signchange_ind(i) == 1            % if we go to pos (so we have been in valley)
        if c_min < max_val_c  && counter > t_crit   % extra discrimination factors (concentration + duration)
            c_val(ind_val) = c_sum/counter;         % average concentration in the valley (mM)
            t_val(ind_val) = t_sum;                 % total time in the valley (s)
        end
        % restart
        c_sum = 0; t_sum = 0; counter = 0; ind_val = ind_val + 1; c_max = min_peak_c; c_min = max_val_c;

    elseif signchange_ind(i) == -1    % we have been in a peak
        if c_max > min_peak_c  && counter > t_crit   % extra discrimination factors (concentration + duration)   
            c_pks(ind_pks) = c_sum/counter;          % average concentration in the valley (mM)   
            t_pks(ind_pks) = t_sum;                  % total time in the valley (s)
        end
        % restart
        c_sum = 0; t_sum = 0; counter = 0; ind_pks = ind_pks + 1; c_max = min_peak_c; c_min = max_val_c;
    end
end
% remove the additional zeros that were allocated, or remain because we were very strict in allocating peaks/valleys)
c_val = nonzeros(c_val); t_val = nonzeros(t_val);
c_pks = nonzeros(c_pks); t_pks = nonzeros(t_pks);
end

%% CSTR ODE model
function [dxdt, q, pode] = fun_CSTR(t,x,p)
% function to be used in ODE solver for calculating concentrations and gas
% mole fractions over time *t*. *x* contains the dissolved compound
% concentrations and gas fractions. *p* indicates whether we are in a 
% peak or valley.

% unpack x
c = x(1:4); y = x(5:7); yh = x(8:10);

D = 1/48;                    % Dilution rate - 1/h
Rrec = 0.76;                  % Biomass recirculation rate
% set characteristics for peak, valley and startup
rpm_su = 75; rpm_low = 20; rpm_high = 950;
vvm = 0.05;

% discriminate between peaks and valleys (matter of seconds)
tsec = t*3600;               % ODE solver time is in hours

% make a pointer pq to assign whether we are in peak, valley or start-up
% phase. This is done by interpolating, to account for transitions between
% peaks and valleys. For proper interpolation, this is only done when 
% tsec > 1 (but we start in a valley anyway). 
if tsec >= 1
    pq = interp1(p,tsec);   % number between 0 (valley) and 1 (peak)
else
    pq = 0;
end
if numel(p) == 2            % in this case p = [0, 0]
    pq = 10;                % arbitrary number
end

% impose the influence of the peak/valleys on kla
if pq < 0.99            % we are in valley!
    [kla_low, PV_low] = fKla(rpm_low, vvm);         % calculate kLa for CO
    rpm = rpm_low;      % set stirrer speed for rest of function
    % make a vector with kla values (which are based on CO)
    kla = [kla_low, kla_low*sqrt(6.01)/sqrt(2.71), kla_low*sqrt(2.56)/sqrt(2.71), 0]; % 1/h [CO, H2, CO2, X]
    pode = 0;           % output for ode solver
    
elseif pq == 10         % start-up phase kla
    [kla_su, PV_su] = fKla(rpm_su, vvm);     rpm = rpm_su;
    kla = [kla_su, kla_su*sqrt(6.01)/sqrt(2.71), kla_su*sqrt(2.56)/sqrt(2.71), 0]; % 1/h [CO, H2, CO2, X]
    pode = 0;           % output for ode solver
    
else                    % we are in peak
    [kla_high, PV_high] = fKla(rpm_high, vvm);   rpm = rpm_high;
    kla = [kla_high, kla_high*sqrt(6.01)/sqrt(2.71), kla_high*sqrt(2.56)/sqrt(2.71), 0]; % 1/h [CO, H2, CO2, X]
    pode = 1;           % output for ode solver
end

% constants for CSTR operation
cin = [0,0,0,Rrec*c(4)];            % Inlet concentrations of dissolved compounds mM [CO, H2, CO2, X]
yin = [0.5, 0.2, 0.3, 0];           % Inlet gas mole fractions [CO, H2, CO2]
T = 273.15 + 37;                    % Temperature - K (37C)
Rgas = 8.314;                       % Universal gas constant - J/K/mol
[H_CO, H_H2, H_CO2] = Henry(T);     % Henry constants (mol/m3/atm)
p_atm = 1;                          % pressure in atm

% Yields for biomass production (Almeida Benalcazar 2020)
Yxco = 0.041;                       % molX/molCO
Yxco2 = 0.021;                      % molX/molCO2
Yxh2 = Yxco2 * 2/6;                 % molX/molH2

% calculate solubility and q-rates
csol = p_atm.*y'.*[H_CO, H_H2, H_CO2]; csol(4) = 0;       % solubility - mM [CO, H2, CO2, X]
q(1) = -1*CO_uptake(c(1));                      % CO uptake rate - molCO/molx/h
q(2) = -1*H2_uptake(c(2), c(1));                % H2 uptake rate - molH2/molx/h
q(3) = q(1)*4/6 - q(2)*1/3;                     % CO2 uptake/production rate - molCO2/molx/h
q(4) = -1*(q(1)*Yxco + q(2)*Yxh2);              % Biomass growth rate - molX/molX/h - 1/h
R = 1*q*c(4);                                   % Rate vector - mol/m3/h

% Retrieve data for gas phase mass balancing (gas inflow [m3/h], gas volume
% in dispersion [m3], gas volume in headspace [m3], volume in liquid [m3],
% total reactor volume [m3] and gas hold-up [m3g/m3D]. 
[Fgin, Vgd, Vgh, Vl, Vr, eg] = Volumes(vvm, rpm);

% liquid phase mass balances (1 for each compound [CO, H2, CO2, X]). 
dcdt = zeros(1,length(c));
for i=1:length(c)
    dcdt(i) = D*cin(i) - D*c(i) + kla(i)*(csol(i)-c(i)) + R(i);         % mol/m3/h
end

% Calculate total amount of gas involved in gas-liquid mass transfer
SMT = kla(1)*(csol(1)-c(1)) + kla(2)*(csol(2)-c(2)) + kla(3)*(csol(3)-c(3));    % mol/m3/h 
p_Pa = 101325;                                  % pressure - Pa
Fgout = Fgin - SMT*Vl*Rgas*T/p_Pa;              % gas ouflow - m3/h

% gas mole balance in dispersion
dydt = zeros(1,3);
for i=1:length(dydt)
    dydt(i) = Fgin/Vgd*(yin(i)-y(i)) + Vl/Vgd*(Rgas*T)/p_Pa*( y(i)*SMT - kla(i)*(csol(i)-c(i)) );
end

% gas mole balance in headspace
dyhdt = zeros(1,3);
for i=1:length(dyhdt)
    dyhdt(i) = Fgout/Vgh*(y(i) - yh(i));
end

% return the ODEs
dxdt = [dcdt, dydt, dyhdt]';

end

%% Peaks development from lifelines
function [t_begin_peak, t_end_peak, valley_times, peak_times] = Peak_development(Npeaks, peak_time_bins_clCO, prob_peak_time_clCO, valley_time_bins_clCO, prob_valley_time_clCO)
% load probability distribution of peaks and valleys from lifelines
% contains the probs and bins for peaks at specific time. 
%
% allocate variables for SD
t_end_peak = 0;                             % times in the peak 
t_begin_peak = zeros(1, Npeaks);            % time points at which we arrive in in peak
valley_times = zeros(1, Npeaks+1);          % vector with times in valleys (one more)
peak_times = zeros(1, Npeaks);              % vector with times in peak
% variables from CFD
p_val = prob_valley_time_clCO;              p_peak = prob_peak_time_clCO;
t_val = valley_time_bins_clCO(1:end-1);     t_peak = peak_time_bins_clCO(1:end-1);

for i = 1:Npeaks
    % choose a "random" value for the valley and peak time using the probability distribution
    while valley_times(i) == 0
        valley_times(i) = randsrc(1,1,[t_val'; p_val' ]);
    end
    while peak_times(i) == 0
        peak_times(i) = randsrc(1,1,[t_peak'; p_peak' ]);
    end
    % update the storage vectors with valley and peak times
    t_begin_peak(i) = t_end_peak(end) + valley_times(i);
    t_end_new = t_begin_peak(i) + peak_times(i);
    t_end_peak(i) = t_end_new;
    
    % the last valley
    if i == Npeaks
        valley_times(i+1) = randsrc(1,1,[t_val'; p_val' ]);
        t_begin_peak(i+1) = t_end_peak(end) + valley_times(i+1);
    end
        
end

end

%% Kullback-Leibler Divergence
function KL_sum = KL_divergence(CFD_probs, SD_data, bins)
    % determine probability SD_data
    h = histogram(SD_data, bins, 'Normalization','Probability', 'Visible', 'off'); % create the histogram for the specific selection of data
    SD_probs = h.Values';            % retrieve the probability for the SD data

    % calculate the KL-divergence as distance from the CFD data
    KLi = CFD_probs .* log((CFD_probs + eps) ./ (SD_probs + eps));
    KL_sum = sum(KLi);              % sum the KL of every bin to return the sum
end

%% Functions to calculate uptake rates, kla, Henry, and volumes

function q = H2_uptake(cH2, cCO)
% inputs are concentrations in mM
qh2max = 2.565;                 % molH2/cmol/h
K_m = 0.025;                    % mol/m3
K_i = 0.025;                    % mol/m3

q = qh2max*cH2/(K_m+cH2)*(1/(1+cCO/K_i)); % mol/molx/h
end


function q = CO_uptake(cCO)
% inputs are concentrations in mM
qcomax = 1.459;                 % molco/cmol/h
K_m = 0.042;                    % mM
K_i = 0.25;                     % mM^2

q = qcomax*cCO/(K_m+cCO+cCO^2/K_i); % mol/molx/h
end


function [H_CO, H_H2, H_CO2] = Henry(Temp)
% Calculates Henry's constants as a fucntion of Temperature - Sander 2015
H_CO = 9.5E-4*exp(1300*(1/Temp - 1/298.15))*1000;   % mol/m3/atm
H_H2 = 7.8E-4*exp(500*(1/Temp - 1/298.15))*1000;    % mol/m3/atm
H_CO2 = 3.4E-2*exp(2400*(1/Temp - 1/298.15))*1000;  % mol/m3/atm

end


function [kla, Pv] = fKla(rpm, vvm)
% reactor characteristics
Vl = 0.002;                 % working volume - m3
Vr = 3e-3;                  % total volume - m3
Dr = 0.13;                  % diamter - m
Ar = 1/4 * pi * Dr^2;       % area - m2
Di = 0.06;                  % impeller diameter - m
Po = 1.5; rhol = 993;       % Power number and density kg-m3
ug = vvm * Vr;              % gas flow - m3/min
Fgin = ug/60;               % gas flow - m3/s
ugs = Fgin/Ar;              % superficial gas flow velocity - m/s

% parameters for kla determination
DLCO = 2.71; DLO2 = 2.8;    % Diffusivities of CO and O2
fkla = 1.5;                 % Broth affects kla with factor fkla
T = 273.15 + 37;            % Temperature - K 

% kla model from Garcia Ochoa and Van 't Riet
Pug = Po*rhol*((rpm/60)^3)*(Di^5);      % ungassed power number
alpha = 0.783; beta = 0.459;            % characteristics for Rusthon turbine
PG = alpha * ( Pug^2*(rpm/60)*Di^3/(Fgin/3600)^(0.56) )^(beta);     % Gassed power number
kla = fkla * (1.022^(T-293.15)) * (2.6e-2*(PG/Vl)^(0.4)*ugs^(0.5)) * sqrt(DLCO/DLO2) * 3600; % kla - 1/h
Pv = PG/Vl;
end


function [Fgin, Vgd, Vgh, Vl, Vr, eg] = Volumes(vvm, rpm)
% reactor characteristics
Vl = 0.002;                 % working volume - m3
Vr = 3e-3;                  % total volume - m3
Dr = 0.13;                  % diamter - m
Ar = 1/4 * pi * Dr^2;       % area - m2
Di = 0.06;                  % impeller diameter - m
ug = vvm * Vr;              % m3/min
Fgin = ug/60;               % m3/s
ugs = Fgin/Ar;              % m/s

% Approximation of gas hold-up using Kudrewiski 1982 or Garcia Ochoa 2004
N = rpm/60;                 % stirrer speed - 1/s
g = 9.81;                   % gravity is always there - m/s2
rho_L = 993;                % density liquid - kg/m3
rho_G = 1.25;               % density gas - kg/m3 
sigma = 0.072;              % surface tension - N/m
A = 0.819 * ugs^(2/3)*N^(2/5)*Di^(4/15)/(g^(1/3)) * (rho_L/sigma)^(1/5) * (rho_L/(rho_L-rho_G)) * (rho_L/rho_G)^(-1/15);
eg = A/(A+1);               % Gas hold up [m3g/m3d]

% retrieve gas volume in dispersion and headspace via volume balances
Vgd = eg*Vl/(1-eg);         % gas in dispersion - m3
Vgh = Vr - Vl - Vgd;        % gas in headspace - m3
Fgin = Fgin * 3600;         % gas inflow - m3/h
end