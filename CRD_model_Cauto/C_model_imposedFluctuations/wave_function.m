
% load data from 150 mM case
load('U:/MicroSynC/Lars/PhD/FLUENT/Kinetic_model/Backflow_bc_update/Ac_Conc_var/No_species_upd/CX_150/Data_lifelines_cx150_500s.mat')
data_150 = [cl_co_s_ll', q_co_s_ll', q_h2_s_ll' q_ac_s_ll', q_etoh_s_ll',c_fd_s_ll', c_fd_std_ll', cl_co_std_ll',cl_h2_s_ll',cl_h2_std_ll',c_donor_s_ll',c_donor_std_ll'];

% get average concentration and its standard deviation
co_conc_av = mean(cl_co_s_ll);
co_conc_std = mean(cl_co_std_ll);
h2_conc_av = mean(cl_h2_s_ll);
h2_conc_std = mean(cl_h2_std_ll);
cfd_conc_av = mean(c_fd_s_ll( find(c_fd_s_ll < 100)) )
cfd_conc_std = std(c_fd_std_ll( find(c_fd_s_ll < 100)) )

oscillation_length = 100;
time_vec = 0:0.1:500;
co_sin = co_conc_std * sin(2*pi/oscillation_length*time_vec) + co_conc_av;
h2_sin = h2_conc_std * sin(2*pi/oscillation_length*time_vec) + h2_conc_av;
plot(time_vec, co_sin, time_vec, h2_sin)