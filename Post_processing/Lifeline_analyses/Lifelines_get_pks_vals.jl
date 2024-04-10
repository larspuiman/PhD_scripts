
# This script analyses the individual lifelines to get peaks and valleys of a certain concentration variable
#   Requires the input of a matrix containing the variable per lifeline and over time


using Statistics
using Random
using DelimitedFiles
using Plots
using BenchmarkTools
using JLD2
using Markdown

"""
    Lifeline_peaks_valleys(dt, vector, min_peak, max_valley )
    
Computes the time and average values during peaks and valleys of a certain (concentration/q) vector from a lifeline observed over time
Requires the input of the time difference between vector values, a vector with values, a minimum value for a peak, and a maximum value to consider being in a valley

""" 
function Lifeline_peaks_valleys(dt, vector, min_peak, max_valley )
 
    # function to analyze lifelines, input is time vector, a vector with variables of interests, concentration/q-rates, minimum values for in peak and maximum values for in the valley
 
    # indexing
    ind_valley, ind_peak = 1,1;                     # indexes for looping
    t_crit = 10;                                    # only store after at least 10 timesteps in peak/valley

    # preallocate memory
    vector_sum, t_sum, counter = 0, 0, 0;
    max_peak, min_valley = min_peak, max_valley;    # to check whether we reached a "real" peak or valley
    t_valley, vec_valley, t_peak, vec_peak = zeros(100,1), zeros(100,1), zeros(100,1), zeros(100,1);

    # 1: check difference with the averages
    vec_av = mean(vector);                          # average value during the lifeline
    vector_averaged = vector .- vec_av;             # vector with average substracted (so that negative values are below average)

    # 2: find where we go from peak to valley (or vice versa, so where is sign-change)
    ind_signchange = sign.(diff(sign.(vector_averaged)))

    # 3: now we are going through the lifeline and assign the data to peaks or valleys
    for vector_index in 1:length(vector)-1

        vector_sum = vector_sum + vector[vector_index];         # make cumulative value of vector (conc, q rate)
        t_sum = t_sum + dt;                                     # make cumulative value of time (duration in peak or valley)
        max_peak = max(max_peak, vector[vector_index]);         # check whether the value in the lifeline is bigger than the minimum value we should have in peak
        min_valley = min(min_valley, vector[vector_index]);     # check whether the value in the lifeline is smaller than the max value we should have in valley
        counter += 1;
        # until here we have cumulative data, which should be corrected for duration. 
      
        
    # 4: we are going to check whether we have been in a peak or valley and save the time and vector values in there.
        if ind_signchange[vector_index] == 1                    # if we go to positive values (so we have been in valley)
            # check whether we have been long enough in valley and whether it was "deep" enough
            if min_valley < max_valley && counter > t_crit
                t_valley[ind_valley] = t_sum;                   # store time in valley
                vec_valley[ind_valley] = vector_sum/counter;    # store vector value in valley
                ind_valley += 1;
            end
            # reset all values for next peak / valley
            vector_sum, t_sum, counter = 0, 0, 0;
            max_peak, min_valley = min_peak, max_valley;    # to check whether we reached a "real" peak or valley

        elseif ind_signchange[vector_index] == -1               # if we have been in peak (then we go negative)
            # check whether we have been long enough in peak abd whether it was "high" enough
            if max_peak > min_peak && counter > t_crit
                t_peak[ind_peak] = t_sum;                   # store time in valley
                vec_peak[ind_peak] = vector_sum/counter;    # store vector value in valley
                ind_peak += 1;
            end
            # reset all values for next peak / valley
            vector_sum, t_sum, counter = 0, 0, 0;
            max_peak, min_valley = min_peak, max_valley;    # to check whether we reached a "real" peak or valley
        end
    end

    return t_valley, vec_valley, t_peak, vec_peak, vec_av;
end

"""
    Lifeline_loader(Lifeline_number)
    Loads the lifeline specified by its number and location (PATH). Lifeline file is an .out file.    
"""
function Lifeline_loader(Lifeline_number::Int64, PATH)
    cd("$PATH")
    Lifeline_data = readdlm("Lifeline_$Lifeline_number.out",',');

    t = Lifeline_data[:,1]
    clCO = Lifeline_data[:,5];
    clH2 = Lifeline_data[:,6];
    clCO2 = Lifeline_data[:,7];
    qCO = Lifeline_data[:,5];
    qH2 = Lifeline_data[:,6];
    qCO2 = Lifeline_data[:,7];
    return t, clCO, clH2, clCO2, qCO, qH2, qCO2
end

"""
    Rolling_average(vector)

    Computes the rolling average of a certain variable over time   
"""
function Rolling_average(vector)
    vec_rav = zeros(length(vector)); vec_rav[1] = vector[1];
    for i in 2:length(vector)
        vec_rav[i] = (vec_rav[i-1]*(i-1) + vector[i])/i
    end
    return vec_rav
end
"""
    Rolling_std(vector, vec_rav)

    Computes the standard deviation of a certain variable over time, requires the vector with the variable as well as its rolling average
"""
function Rolling_std(vector, vec_rav)
    vec_rvar = zeros(length(vector));
    vec_rstd = zeros(length(vector));
    M2 = zeros(length(vector));
    for i in 2:length(vector)
        M2[i] = M2[i-1] + ((vector[i] - vec_rav[i-1])*(vector[i] - vec_rav[i]));  # rolling variance
        vec_rvar[i] = M2[i]/i;
    end
    vec_rstd = sqrt.(vec_rvar)

    return vec_rstd
end


# VARIABLE TO BE READ AND ITS LOCATION
variable = "clCO"; cx = 10;
if cx == Int(10); run_ind = 5; end;
if cx == Int(25); run_ind = 3; end;
if cx == Int(5); run_ind = 1; end;
cd("U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run$run_ind-DynField_$cx-gl/Lifelines")

# HOW MANY LIFELINES SHOULD BE READ? AND FOR HOW LONG? Disregard the first mixing time?
max_lifelines = 41715;                                     # amount of lifelines (can be more than # in file)
tmix = 90;                                                  # mixing time of reactor (s)
t_start = 3000; dt = 0.1; t_end = 3650;                     # start-time of DPM method, time step, and end-time.
t = t_start:dt:t_end;                                       # time-vector for the DPM method                                    
t_total = Vector((t_start+tmix):dt:t_end);                  # time-vector that should be saved (does not include the 1st mixing time)
t_total_ind = Vector(tmix/dt:1:(length(t)-1));              # indices of stored-time vector

# selective time data for storing rolling average and standard deviations
time_back = 50; freq = 1;                                        # time_back: how often should data be stored? every 50s. Freq: if this time should be higher or lower                 
tm_back_rel = Vector(0:freq:(length(t_total)/(time_back)*dt));   # amount of datapoints and their locations stored
tm_back_abs = (tm_back_rel .* time_back) .+ (t_start + tmix);    # time value of these data points
tm_back_abs_ind = zeros(length(tm_back_abs));                    # retrieve the indices for these datapoints
for i in 1:length(tm_back_abs)
    a = findall(x -> x == tm_back_abs[i], t_total)
    tm_back_abs_ind[i] = a[1];
end
tm_back_abs_ind = floor.(Int, tm_back_abs_ind);                   # indices for all datapoints that should be stored 

# EULERIAN MEAN VALUES (tecplot @ 3000 s (200 s time average))
clCO_Emean = [6.99E-01, 7.82E-02, 3.10E-02, 2.02E-02, 4.74E-03];              # mM
clH2_Emean = [2.93E-01, 3.69E-02, 0.01154882, 6.84E-03, 9.27E-04];            # mM
qCO_Emean = [3.92E-01, 6.13E-01, 4.52E-01, 3.57E-01, 1.33E-01];               # mol/mol/h
qH2_Emean = [8.82E-02, 2.62E-01, 2.02E-01, 1.59E-01, 5.99E-02];               # mol/mol/h
ind_Emean = 0;
if cx == 2.5; ind_Emean = 1; end;
if cx == 5.0; ind_Emean = 2; end;
if cx == 7.5; ind_Emean = 3; end;
if cx == 10; ind_Emean = 4; end;
if cx == 25; ind_Emean = 5; end;
E_clCO_av = clCO_Emean[ind_Emean];                                          # mmol/m3 or mM
E_clH2_av = clH2_Emean[ind_Emean];                                          # mmol/m3 or mM
E_qCO_av = qCO_Emean[ind_Emean];                                            # mol/mol/h
E_qH2_av = qH2_Emean[ind_Emean];                                            # mol/mol/h

# WRITE THE RESULTS IN A jld2 DATAFILE?
write_data = 1;

# set general description for min and max value in peak based upon Eulerian average
factor_pkval = 1.5;                                 # 1.5 for 10 and 25 g/L and 2 for 5 g/L
if variable == "clCO"
    min_peak = E_clCO_av*factor_pkval;
    max_valley = E_clCO_av/factor_pkval; 
elseif variable == "clH2"
    min_peak = E_clH2_av*factor_pkval;
    max_valley = E_clH2_av/factor_pkval; 
elseif variable == "qCO"
    min_peak = E_qCO_av*factor_pkval;
    max_valley = E_qCO_av/factor_pkval; 
elseif variable == "qH2"
    min_peak = E_qH2_av*factor_pkval;
    max_valley = E_qH2_av/factor_pkval; 
end


@time begin
    # allocate memory for the lifeline analysis results
    t_valley, vec_valley, t_peak, vec_peak, vec_av = zeros(100, max_lifelines), zeros(100, max_lifelines), zeros(100, max_lifelines), zeros(100, max_lifelines), zeros(max_lifelines);
    vec_rtav = zeros(length(tm_back_abs), max_lifelines); vec_rtstd = zeros(length(tm_back_abs), max_lifelines); vec_rtcov = zeros(length(tm_back_abs), max_lifelines);
    global iter = 1; 
    global iter_aborted = 0; 

    # open the .out file with the variable data
    fp = open("sorted_lifelines/$variable.out") 

    # read file line-by-line
    for (i, line) in enumerate(eachline(fp))
        # as long as we are not at the maximum amount of lifelines that should be plotted
        if iter <= max_lifelines
            # obtaind data from line by splitting it and converting each value in the string to a Float 
            a = split(line, "\t")
            array = [parse(Float64,value) for value in a ]
            
            # The first value is the particle number, the later values are the time data (10000 points)< vector with q or c values
            vector = array[(floor(Int,t_total_ind[1]+1)):end]; 

            # if lifeline is aborted during simulation (then the length of vector is too short), so disregard the lifeline
            if length(vector) < length(t_total)-100
                global iter_aborted += 1;                                                      # to store how many lifelines are aborted     
                println("Iter = $iter, Aborted lifeline!");        
                # go to the next lifeline
                global iter += 1;      
                continue        
            end

            # now retrieve the peaks and valleys of our lifeline and store them for this specific lifeline.
            t_valley[:,iter], vec_valley[:,iter], t_peak[:,iter], vec_peak[:,iter], vec_av[iter] = Lifeline_peaks_valleys(dt, vector, min_peak, max_valley); 

            # determine rolling average and standard deviation
            vec_rtav_temp = Rolling_average(vector);
            vec_rtstd_temp = Rolling_std(vector, vec_rtav_temp);
            vec_rtcov_temp = vec_rtstd_temp./vec_rtav_temp;
            # only store rolling values at time points we want
            vec_rtav[:,iter] = vec_rtav_temp[tm_back_abs_ind];
            vec_rtstd[:,iter] = vec_rtstd_temp[tm_back_abs_ind];
            vec_rtcov[:,iter] = vec_rtcov_temp[tm_back_abs_ind];


            if iter%1000 == 0
                println("Lifeline $iter analysed")
            end
            # go to the next lifeline
            global iter += 1;

        else
            # now we reached the amount of lifelines we wanted, so can we stop the loop.
            break
        end
    end
    print("escaped!")    
end #time

# # put data into arrays nd remove the zeros
t_valley = t_valley[:]; t_peak = t_peak[:]; vec_valley = vec_valley[:]; vec_peak = vec_peak[:];
deleteat!(t_peak, t_peak .== 0); deleteat!(vec_peak, vec_peak .== 0);
# find zeros for the valleys. Delete when time in valley = 0 (as there could be valleys with 0 conc but not with 0 time)
ind_valley_del = (t_valley .== 0);
deleteat!(t_valley, ind_valley_del .== 1); deleteat!(vec_valley, ind_valley_del .== 1);
# deleteat!(vec_av, vec_av .== 0);


if write_data == 1
    jldsave("sorted_lifelines/Peaks_valleys_av_$variable-$cx.jld2"; iter, t_valley, vec_valley, t_peak, vec_peak, vec_av);
    jldsave("sorted_lifelines/Rolling_tav_std_$variable-$cx.jld2"; iter, vec_rtav, vec_rtstd);
end


# vec_rtn_av = zeros(size(vec_rtav)); vec_rtn_std = zeros(size(vec_rtstd)); vec_rtn_cov = zeros(size(vec_rtcov));
# for i = 1:size(vec_rtav,1)
#     vec_rtn_av[i,:] = Rolling_average(vec_rtav[i,:]);
#     vec_rtn_std[i,:] = Rolling_std(vec_rtav[i,:], vec_rtn_av[i,:]);
#     vec_rtn_cov[i,:] = vec_rtn_std[i,:]./vec_rtn_av[i,:];
# end

# # sample variance 
# lang_average = mean(vec_av);

