clearvars

conc_act = load('Ini_conc.dat');
conc_ref = load('Ref_conc.dat');

t0 = 0; tend = 50/3600; dt = 1e-4;
num_tsteps = (tend)/dt;
tspan = [t0:0.1/3600:tend];
tstart = tic;
opt = odeset('AbsTol',1e-3,'RelTol',1e-3);
[t, c] = ode15s(@dcdt_fun,tspan,conc_act, opt);
t = t*3600;



%%
save('Batch_matlab.mat','t', 'c');





%% plot data using subplot in for loop
figure()
compounds = {'CO', 'H_2','CO_2', 'AcT', 'EtOH', 'BDO', 'AcCoA', 'NADH', ...
    'NADPH', '[Ac^-]_{IC}', 'For^{-}', 'Fd_{red}^{2-}'};

for i = 1:size(c,2)
    subplot(4,3,i)
    plot(t, c(:,i), 'r-', 'LineWidth',2); hold on;
    ylabel('{\itc_{i}} (mol m^{-3})');
    
    title(compounds{i});
    xlabel('time (s)'); xlim([0, 20]);

end

set(gcf,'position', [117.8 122.6 1216 640.4]);

%% compare dcdts
[dc_dt, dc_dt_original, rates_new] = dcdt_fun(t, c);
figure()
for i = 1:size(c,2)
    subplot(4,3,i)
    plot(t, dc_dt(:,i), 'r-', 'LineWidth',2); hold on;
    plot(t, dc_dt_original(:,i), 'b--', 'LineWidth',1.5); hold on;
    title(compounds{i});
    xlabel('time (s)'); xlim([0, 20]);
    ylabel('{\it dc_{i}/dt} (mol m^{-3} h^{-1})')
    
   
end
    
set(gcf,'position', [117.8 122.6 1216 640.4]);
  %% CALCULATE RATES
 
fluxes = {'FDH', 'ACAS', 'CODH', 'HYD', 'AcS', 'EtS', 'BDOS','Nfn','Rnf','Ana','AcX', 'Import'};

% 
% %%
% figure()
% for i = 1:12
%     subplot(4,3,i)
% %     plot(t, rates_cfd(:,i)); hold on;
%     plot(t_data_F, rates_conc(:,i)); hold on;
%     plot(t_data_F, rates_new(:,i));
%     std_rate = std(rates_conc(:,i));
%     av_rate = mean(rates_conc(:,i));
% %     plot(t, rates_conc_constr(:,i));
% %     ylim([av_rate-3*std_rate, av_rate+3*std_rate]); 
%     title(fluxes{i});
%     xlabel('time (s)');     
%     ylabel('mol/molx/h')
% end
%     set(gcf,'position', [117.8 122.6 1216 640.4]);








    

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




%%
%% calculate dcdt
function [dc_dt_update, dc_dt, rates_new, imax] = dcdt_fun(t, conc_act)

% if we only have 1 moment in time, then transpose the DATA
if size(conc_act,2) == 1
    conc_act = conc_act';
end


% get rates
rates = rates_fun_v(conc_act);
CX = 200; Vx = 5.89e-5;                         % mol/m3 liq, and m3 liq / mole cell
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
    
% if we only have 1 moment in time, then transpose the DATA
if size(dc_dt,1) == 1
    dc_dt = dc_dt';
    dc_dt_update = dc_dt_update';
end
    
    
end
