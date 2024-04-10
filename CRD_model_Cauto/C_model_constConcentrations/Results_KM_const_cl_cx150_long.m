clearvars
format compact
close all


conc_array =  [ 0.05:0.01:0.14];
q_rates_mean = zeros(length(conc_array), 6);
conc_end = zeros(length(conc_array), 12);

bal = zeros(length(conc_array), 1);
Ac_2_EtOH = zeros(length(conc_array), 1);
CO_2_EtOH = zeros(length(conc_array), 1);
qCO_mean = zeros(length(conc_array), 1);
for i = 1:length(conc_array)
    conc_Ac = conc_array(i);
    [conc_end(i,:), q_rates_mean(i,:) ] = load_data(conc_Ac);
end

%%
figure()
plot(conc_array, q_rates_mean(:,1))
xlabel('c AcT')
ylabel('qCO (mol/mol/h)')

%% check for balance

figure()
subplot(211)
plot(conc_array, -q_rates_mean(:,5)./q_rates_mean(:,1), 'b-', 'LineWidth',2) ; hold on;
xlabel('c_{CO} (mM)')
ylabel('Y_{EtOH/CO} (mol_{EtOH}/mol_{CO})')

subplot(212)
plot(conc_array, (q_rates_mean(:,4)./q_rates_mean(:,5)), 'b-', 'LineWidth',2) ; hold on;
xlabel('c_{CO} (mM)')
ylabel('Y_{Ac/EtOH} (mol_{Ac}/mol_{EtOH})')

%% save data
save('Results_cx150mM_long/Results_KM_const_cl_cx150.mat', 'conc_end', 'q_rates_mean')



function [conc_end, q_rates_mean ] = load_data(conc_CO)

    conc_data_str = ['Results_cx150mM_long/conc_data_cCO_' num2str(conc_CO) '.dat'];
    dcdt_data_str = ['Results_cx150mM_long/dcdt_data_cCO_' num2str(conc_CO) '.dat'];
    t_data_str = ['Results_cx150mM_long/t_data_cCO_' num2str(conc_CO) '.dat'];

    conc_data = load(conc_data_str)';
    dcdt_data = load(dcdt_data_str)';
    t_data = load(t_data_str);

    tmax = 1000;
    tmin = 100;
    ind_max = find(t_data > tmin & t_data < tmax);
    conc_data = conc_data(ind_max, :);
    dcdt_data = dcdt_data(ind_max, :);
    t_data = t_data(ind_max, :);
%     
    
    compounds = {'CO', 'H2','CO2', 'AcT', 'EtOH', 'BDO', 'AcCoA', 'NADH', ...
        'NADPH', '[Ac-]i', 'For-', 'Fdred2-'};

    %% calculate dcdt
    [dcdt, dcdt_update, rates_new, imax] = dcdt_fun(conc_data);
    figure(1)
    for i = 1:size(conc_data,2)
        subplot(4,3,i)
        plot(t_data, conc_data(:,i)); hold on;
        ylabel('c (mM)')
    %     yyaxis right
    %     plot(t_data, dcdt_data(:,i));
    %     ylabel('dcdt (mol/m3/h)')
        title(compounds{i});
        xlabel('time (s)')
        set(gcf,'position', [117.8 122.6 1216 640.4]);
    end




    %saveas(gcf,'Results/conc_devs_tlang_1ms_1000molm3.png')
      %% CALCULATE RATES

    fluxes = {'FDH', 'ACAS', 'CODH', 'HYD', 'AcS', 'EtS', 'BDOS','Nfn','Rnf','Ana','AcX', 'Import'};

    rates_conc = zeros(size(conc_data,1), 12);
    [rates_conc] = rates_fun_v(conc_data);
  
%     figure(3)
%     for i = 1:size(conc_data,2)
%         subplot(4,3,i)
%         plot(t_data, rates_new(:,i)); hold on;
%         ylabel('rate (mol/mol/h)')
%     %     yyaxis right
%     %     plot(t_data, dcdt_data(:,i));
%     %     ylabel('dcdt (mol/m3/h)')
%         title(fluxes{i});
%         xlabel('time (s)')
%         set(gcf,'position', [117.8 122.6 1216 640.4]);
%     end

    %% calculate balance
    CX = 150; Vx = 5.89e-5;          
    qCO = dcdt_data(:,1)./CX;
    qH2 = dcdt_data(:,2)./CX;
    qCO2 = dcdt_data(:,3)./CX;
    qAc = dcdt_data(:,4)./CX;
    qEtOH = dcdt_data(:,5)./CX;
    qBDO = dcdt_data(:,6)./CX;
    mu = rates_new(:,10);
    qAcCoA_in = dcdt_data(:,7).*Vx;
    qAc_in = dcdt_data(:,10).*Vx;
    qFor_in = dcdt_data(:,11).*Vx;
    qNADH = dcdt_data(:,8).*Vx;
    qNADPH = dcdt_data(:,9).*Vx;
    qFD = dcdt_data(:,12).*Vx;
    av_rates = mean(rates_new);
    
    % 
    c_rates = [qCO, qCO2, qAc*2, qEtOH*2, qBDO*4, mu, qAcCoA_in*2, qAc_in*2, qFor_in];
    c_rates_mean = mean(c_rates);
    c_rates_mean = [c_rates_mean(1:6), sum(c_rates_mean(7:end))];
    bal = sum(c_rates_mean);
    bal_rel = bal/(abs(c_rates_mean(1) + c_rates_mean(3)));

    Ac_2_EtOH = mean( (abs(qAc) )./abs(qEtOH));
    CO_2_EtOH = mean(abs(qCO)./abs(qEtOH));
%     CO_2_EtOH = (mean(abs(qCO))./mean(abs(qEtOH)));
    qCO_mean = mean(abs(qCO));
    q_rates_mean = mean([qCO, qH2, qCO2, qAc, qEtOH, mu]);
    conc_end = conc_data(end,:);
    q_rates_end = [qCO(end), qH2(end), qCO2(end), qAc(end), qEtOH(end), mu(end)];


    
    e_rates = [qCO*2, qH2*2, qAc*8, qEtOH*12, qBDO*22, mu*4.2, qAcCoA_in*9, qAc_in*8, qFor_in*2, qNADH*2, qNADPH*2, qFD*2];
    e_rates_mean = mean(e_rates);
    e_rates_mean = [e_rates_mean(1:6), sum(e_rates_mean(7:end))];

end



%%
function rates = rates_fun_v(conc_act)
    conc_act = real(conc_act);
    % constants
    pKa_HAc = 4.756; pH_EC = 5;
    cell_surf_area = 321;                     % surface area of cells (m2 / molx) */
    k_HAcD = 3.85e-5;                         % diffusivity of acetic acid (m / h) */

    % load data
    elasticity_matrix = load('eps_matrix.dat');
    rates_IC_ref = load('Ref_rates.dat');
    conc_ref = load('Ref_conc.dat');
    min_concv = load("Min_conc.dat");

    % calculate derivatives
    neg_conc = zeros(1,length(conc_ref));
    % check whether the concentrations are below 0 and store when that is
    rates_llcor = zeros(size(conc_act,1), length(rates_IC_ref));
    
    % check whether the concentrations are below 0 and store when that is

    
for k = 1:size(conc_act,1) % temporal domain 
    for i = 1:length(rates_IC_ref)
        sumelast = 0; 
        for j = 1:length(conc_ref)
            sign = 1;
             if conc_act(k,j) < min_concv(j)
                 conc_eff = min_concv(j);
             else
                 conc_eff = conc_act(k,j);
             end
%              conc_eff = conc_act(k,j);
%             if (j == 12 && conc_act(k,j) < 0.01)
%                 conc_eff = 0.01;
%             end
            logterm = log(conc_eff/conc_ref(j));
%             if ((elasticity_matrix(i,j) < 0) && (logterm < -2) && (i ~= 2)  ) 
%                 sign = -1;
%             end

            sumelast = sumelast + elasticity_matrix(i,j)*logterm*sign;
        end
        rates_llcor(k,i) = 1 + sumelast;
        if (rates_llcor(k,i) < 0 && ( i == 2 || i == 5 || i == 9 || i == 10 || i == 11 ))
            rates_llcor(k,i) = 0;
        end
    end
end
    rates_llcor = rates_llcor';
    rates = rates_IC_ref.*rates_llcor;         % (mol/mol/h)
    % calculate rate of acetate back-diffusion (mol/mol/h)
    c_HAc_EC = (conc_act(:,4)*10^pKa_HAc)/( 10^(pH_EC)+ 10^(pKa_HAc));    % c HAc out (mol/m3)
    r_HAcD = cell_surf_area*k_HAcD*c_HAc_EC;  % back diffusion rate (mol/molx/h) */

    rates(12,:) = r_HAcD;
    rates = rates';
end


%% calculate dcdt
function [dc_dt, dc_dt_update, rates_new, imax] = dcdt_fun(conc_act)
% get rates
rates = rates_fun_v(conc_act);
CX = 150; Vx = 5.89e-5;                         % mol/m3 liq, and m3 liq / mole cell
    pKa_HAc = 4.756; pH_EC = 5;
    cell_surf_area = 321;                     % surface area of cells (m2 / molx) */
    k_HAcD = 3.85e-5;                         % diffusivity of acetic acid (m / h) */


% calculate rates
mu = rates(:,10); r_HAcD = rates(:,12);
dc_dt(:,1) = -1*(rates(:,2) + rates(:,3))   * CX ;              % CO uptake rate (mol/m3/h) (negative for consumption) */
dc_dt(:,2) = -2*(rates(:,4))                 * CX ;              % H2 uptake rate (mol/m3/h) (negative for consumption) */
dc_dt(:,3) =    (rates(:,3) - rates(:,1))   * CX ;              % CO2 uptake rate (mol/m3/h) (negative for consumption)*/
dc_dt(:,4) =    (rates(:,11) - r_HAcD)       * CX ;              % Net acetate excretion rate (mol/m3/h) */
dc_dt(:,5) =    (rates(:,6))                 * CX ;              % EtOH prod rate (mol/m3/h) */
dc_dt(:,6) =    (rates(:,7))                 * CX ;              % BDO prod rate (mol/m3/h) */
% only compounds that remain in cell (mol/m3_cell/h) */
dc_dt(:,7) =  (rates(:,2) -   rates(:,5) -    2*rates(:,7) - 1/2*rates(:,10)) ...
    .* 1/Vx - mu.*conc_act(:,7);   % AcCoA balance (mol/m3c/h) */
dc_dt(:,8) =  (rates(:,9) - 2*rates(:,2) -     rates(:,6)  - 1/2*rates(:,7) - rates(:,8))...
    .* 1/Vx - mu.*conc_act(8);   % NADH balance (mol/m3c/h) */
dc_dt(:,9) =  (rates(:,4) + 2*rates(:,8) - 1/2*rates(:,1)  -     rates(:,2)) ...
    .* 1/Vx - mu.*conc_act(9);   % NADPH balance (mol/m3c/h) */
dc_dt(:,10) =  (rates(:,5) +   r_HAcD   -     rates(:,6)  -     rates(:,11)) ...
    .* 1/Vx - mu.*conc_act(:,10);   % Ac_i_in balance (mol/m3c/h) */
dc_dt(:,11) = (rates(:,1) -   rates(:,2))  ...
    .* 1/Vx - mu.*conc_act(:,11);  % Formate balance (mol/m3c/h) */
dc_dt(:,12) = (-1/2*rates(:,1) + rates(:,2) + rates(:,3) + rates(:,4) - ...
        rates(:,6) - rates(:,7) - rates(:,8) - rates(:,9)) ...
        .* 1/Vx - mu.*conc_act(:,12);

% min_concv =  [1e-10, 1e-10, 0.1, 0.1, 0.1, 0.1, 1e-9, 1e-12, 1e-12, 1e-12, 1e-9, 1e-15];
min_concv = load("Min_conc.dat")';

max_rate(:,1:6) = (conc_act(:,1:6) - min_concv(1:6)*0.9)./( 1e-4/3600 );
max_rate(:,7:12) = (conc_act(:,7:12) - min_concv(7:12)*0.9)/( 1e-4/3600 );
% max_rate = (conc_act - min_concv*0.9)/( 1e-4/3600 );

[max_rate_i, imax] = max((dc_dt - max_rate)');
dc_dt_update = dc_dt;
for i = 1:size(conc_act,1)
    
    if any(dc_dt(i,:) > max_rate(i,:))
        ref_dcdt = max_rate(i, imax(i));

        dc_dt_update(i,:) = ref_dcdt/ abs( dc_dt(i,imax(i)) ) * dc_dt(i,:);
    end
end

%% recalculate rates from update
A = [0	-1	-1	0	0	0	0	0	0	0	0           0 ;
    0	0	0	-2	0	0	0	0	0	0	0           0 ;
    -1	0	1	0	0	0	0	0	0	0	0           0 ; 
    0	0	0	0	0	0	0	0	0	0	1           -1 ;
    0	0	0	0	0	1	0	0	0	0	0           0 ;
    0	0	0	0	0	0	1	0	0	0	0           0 ; 
    0	1	0	0	-1	0	-2	0	0	-0.5	0       0 ;
    0	-2	0	0	0	-1	-0.5	-1	1	0	0       0 ;
    -0.5	-1	0	1	0	0	0	2	0	0	0       0 ;
    0	0	0	0	1	-1	0	0	0	0	-1          1;
    1	-1	0	0	0	0	0	0	0	0	0           0;
    -0.5	1	1	1	0	-1	-1	-1	-1	0	0       0 ];
A = A(:,1:11);
A(1:6,:) = A(1:6,:).*CX;
A(7:12,:) = A(7:12,:)./Vx;
rates_new = zeros(size(conc_act,1), 12);

for i = 1:size(conc_act,1)
    At = A;
    At(7:12,10) = At(7:12,10) - conc_act(i,7:12)';
    rates_new(i,1:11) = At\dc_dt_update(i,:)';
    % calculate rate of acetate back-diffusion (mol/mol/h)
    c_HAc_EC = (conc_act(i,4)*10^pKa_HAc)/( 10^(pH_EC)+ 10^(pKa_HAc));    % c HAc out (mol/m3)
    r_HAcD = cell_surf_area*k_HAcD*c_HAc_EC;  % back diffusion rate (mol/molx/h) */
    rates_new(i,12) = r_HAcD;
    
    rates_new(i,11) = r_HAcD + rates_new(i,11);
end
    
    
    
end