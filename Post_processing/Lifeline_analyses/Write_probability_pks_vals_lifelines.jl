using JLD2
using StatsBase
using LinearAlgebra
using MAT

function read_data(cx, var)
    if cx == Int(10); run_ind = 5; end;
    if cx == Int(25); run_ind = 3; end;
    if cx == Int(5); run_ind = 1; end;
    # WORKING DIRECTORY OF DATASET
    #cd("U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run$run_ind-DynField_$cx-gl/Lifelines/sorted_lifelines")
    cd("C:/Users/larspuiman/Documents/PhD/BC_Gaslift_modelling/paper_scale_down")
    f = jldopen("data/Peaks_valleys_av_$var-$cx.jld2");
    t_peak = f["t_peak"]
    c_peak = f["vec_peak"]
    t_valley = f["t_valley"]
    c_valley = f["vec_valley"]
    N_lifelines = f["iter"];                            # Number of lifelines from input file (can also be approximation)

    println("$var data is read!")
    return t_peak, c_peak, t_valley, c_valley, N_lifelines;
end

function get_probs(vector)
    max_bin_value = ceil(maximum(vector), sigdigits=2);     # maximum value 
    if max_bin_value > 150; max_bin_value = 150; end;
    min_bin_value = floor(minimum(vector), sigdigits=2);    # minimum value 
    binsize = (max_bin_value-min_bin_value)/100;            # 100 values
    bins = min_bin_value:binsize:max_bin_value;             # the number of bins

    # retrieve probability per bin of the whole dataset
    h1 = fit(Histogram, vector, bins)
    h = normalize(h1, mode=:probability)
    probs = h.weights;
    return probs, bins
end

# SELECT DATA TO READ
cx = 5; var = "clH2";
t_peak, c_peak, t_valley, c_valley, N_lifelines = read_data(cx, var);

# get probabilities for peaks 
p_t_peaks, t_bins_pk = get_probs(t_peak);
p_t_valley, t_bins_val = get_probs(t_valley);
p_c_peaks, c_bins_pk = get_probs(c_peak);
p_c_valley, c_bins_val = get_probs(c_valley);

# file = matopen("prob-peaks-val-$var-$cx-gl.mat", "w")
cd("C:/Users/larspuiman/Documents/PhD/BC_Gaslift_modelling/paper_scale_down/code/scale_down/Results_matlab")
matwrite("prob-peaks-val-$var-$cx-gl.mat", Dict(
            "prob_peak_time_$var" => p_t_peaks,
            "peak_time_bins_$var" => collect(t_bins_pk), 
            "prob_valley_time_$var" => p_t_valley,
            "valley_time_bins_$var" => collect(t_bins_val), 
            "prob_peak_conc_$var" => p_c_peaks,
            "peak_conc_bins_$var" => collect(c_bins_pk), 
            "prob_valley_conc_$var" => p_c_valley,
            "valley_conc_bins_$var" => collect(c_bins_val) ); compress = true)
#close(file)
