using JLD2
using StatsBase
using LinearAlgebra
using PlotlyJS
using Plotly
using DataFrames
using LaTeXStrings

# Plotly.signin("larspuiman", "n9PrIux9FAb530P2aDd6")

# DATASET TO BE ANALYZED
var = "clCO"; cx = 5;
read_data = 1
if read_data == 1
    # WORKING DIRECTORY OF DATASET
    # cd("U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run5-DynField_$cx-gl/Lifelines/sorted_lifelines")
    cd("U:/MicroSynC/Lars/Data_models_thesis/Ch4_CFD_SD_paper/data")   # working directory
    f = jldopen("Peaks_valleys_av_$var-$cx.jld2");
    t_peak = f["t_peak"]
    c_peak = f["vec_peak"]
    t_valley = f["t_valley"]
    c_valley = f["vec_valley"]

    println("$var data has been read!")
end

# STORING YOUR RESULTS
vec_string_plot = "peak";                           # for the plot
letter = "a"; compound = "CO"
save_figs = 1;                                      # switch to save the figures or not
cd("U:/MicroSynC/Lars/Data_models_thesis/Ch4_CFD_SD_paper")   # working directory


function f_probabilities(time_val, time_bins, vec_val, vec_bins)
    # determine probability for the time data
    time_hist = fit(Histogram, time_val, time_bins);            # create the histogram for the time data
    time_hist_p = normalize(time_hist, mode=:probability);      # normalize the histogram based on probability (Sum all p = 1)
    time_prob = time_hist_p.weights;                            # retrieve the probability of the current sample 

    # determine probability for the vector data (Concentration, q-rate)
    vec_hist = fit(Histogram, vec_val, vec_bins);               # create the histogram for the vector data
    vec_hist_p = normalize(vec_hist, mode=:probability);        # normalize the histogram based on probability (Sum all p = 1)
    vec_prob = vec_hist_p.weights;                              # retrieve the probability of the current sample 
 
    # determine probabilities for the combination of time and vector values for each peak or valley. 
    #   we loop through the whole vector of values (all peaks/valleys)
    #   then we check in which bin the current peak/valley is and retrieve its probability. 
    #   after that, we multiply the value of the individual peak and valley bins.
    #       for example, the first peak his time value is in the 20th time bin. 
    #       then we need to store the probability of the 20th time-bin in our temporary storage vector 
    #       but the Concentration value of the same pake is in the 5th concentration bin, so we need to store that value
    #       after that we multiply the individual time + conc probabilities for each peak/valley

    # allocate memory for all peak/valley values
    time_prob_val, vec_prob_val = zeros(length(time_val)), zeros(length(vec_val));
    dtbins = time_bins[2]-time_bins[1]; dvbins = vec_bins[2]-vec_bins[1];

    for i = 1:length(time_val)
        for j = 1:length(time_prob)                 
            # if our time value is in the bin, we need to save the probability that applies to this peak/valley
            if (time_val[i] > time_bins[j]) && (time_val[i] < (time_bins[j] + dtbins))
                time_prob_val[i] = time_prob[j];
            end
        end
        for j = 1:length(vec_prob)
            # if our vector value is in the bin, we need to save the probability that applies to this peak/valley
            if (vec_val[i] > vec_bins[j]) && (vec_val[i] < (vec_bins[j] + dvbins))
                vec_prob_val[i] = vec_prob[j];
            end
        end
    end

    # get the individual probabilities and normalize them based on sum p = 1 
    prob_val_t = time_prob_val .* vec_prob_val;
    prob_val = prob_val_t./sum(prob_val_t);
    return prob_val, time_prob_val, vec_prob_val;
end



# SET PARAMETERS FOR LIFELINE ANALYSIS
N_lifelines = f["iter"];                                        # Number of lifelines from input file (can also be approximation)
N_peaks, N_valleys = length(t_peak), length(t_valley);          # Number of peaks and valleys in the analysis
N_analyzed = Int(5e5);                                           # How many peaks/valleys do we want to use in the analysis?

# SELECT TIME DATA TO BE ANALYZED    
time_values = t_peak[1:N_analyzed];                           # !! vector with the variable that should be analyzed
binsize = 1;                                                    # the width of a bin (1s, 1mM)?
max_bin_value_t = ceil(maximum(time_values), sigdigits=1);      # maximum value (800 s / 80 s for valley / peak)
time_bins = 1.0:binsize:max_bin_value_t;                        # the number of bins

# SELECT VALUE DATA TO BE ANALYZED (CONC OR Q-RATE)
vec_values = c_peak[1:N_analyzed];                            # !! vector with the variable that should be analyzed
min_bin_value = floor(minimum(vec_values), sigdigits = 1);      # minimum value (0.08 mM / 0.5 mM for valley / peak)
max_bin_value_v = ceil(maximum(vec_values), sigdigits = 1);     # maximum value (0.08 mM / 0.5 mM for valley / peak)
binsize = max_bin_value_v/100;                                    # the width of a bin, 100 values
if min_bin_value < binsize; min_bin_value = 0; end              # if the minimum value is super low, then we can beter take 0.
vec_bins = min_bin_value:binsize:max_bin_value_v;               # the number of bins


# retrieve probabilities using the written function
prob_val, time_prob_val, vec_prob_val = f_probabilities(time_values, time_bins, vec_values, vec_bins);

# plot probabilities

# set range of x axis
x_lim = 50; if vec_string_plot == "valley"; x_lim = 200; end

# make layout
layout = Layout(
    title=attr(text="<b>$letter)</b>", x = 0.01, y = 0.9), 
    xaxis_title="<i>t</i><sub>$vec_string_plot</sub> (s)",
    xaxis=attr(range=[0, x_lim], showgrid=true, gridwidth=0.3, gridcolor="lightgrey", 
        zeroline=true, zerolinewidth=1, zerolinecolor="black",
        ticks="inside", tickwidth=2, tickcolor="black", ticklen=5),
    yaxis_title="<i>c<sub>L,$compound,</i></sub><sub>$vec_string_plot</sub> (mol m<sup>-3</sup>)",
    yaxis=attr(range=[0, max_bin_value_v], showgrid=true, gridwidth=0.3, gridcolor="lightgrey",
        zeroline=true, zerolinewidth=1, zerolinecolor="black", 
        ticks="inside", tickwidth=2, tickcolor="black", ticklen=5, 
        automargin=true),
    legend_title="Legend Title",
    font=attr(
        family="times",
        size=26,
        color="black"),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    width = 700, height = 500 );

s = plot(scattergl(
    x = time_values,
    y = vec_values, 
    mode="markers",
    marker=attr(size=2, color=prob_val, showscale=true, line_width=0, colorbar_title = "<i> p </i>", 
    colorscale="Viridis", colorbar_exponentformat="power", cmin=0, cmax = 4e-6)
), layout)

#Plotly.post(s, filename="scatter_$var-$vec_string_plot", auto_open = true)

  
if save_figs == 1
    PlotlyJS.savefig(s,"results/scatter_prob/Figure4$letter-$compound-$vec_string_plot-$cx.png", scale=4);
    PlotlyJS.savefig(s,"results/scatter_prob/Figure4$letter-$compound-$vec_string_plot-$cx.pdf");
    PlotlyJS.savefig(s,"results/scatter_prob/Figure4$letter-$compound-$vec_string_plot-$cx.html");
 #   PlotlyJS.savefig(s,"results/scatter_prob/scatter_$var-$vec_string_plot.eps");
end
display(s)

# s1 = plot(scatter(
#     x = time_values,
#     y = vec_values, 
#     mode="markers",
#     marker=attr(size=2, color=prob_val, showscale=true, line_width=0, colorbar_title = "<i> p </i>", 
#     colorscale="Viridis", colorbar_exponentformat="power")
# ), layout)
# close(s1)
# if save_figs == 1
#     PlotlyJS.savefig(s1,"results/scatter_prob/scatter_$var-$vec_string_plot-s1.png");
#     PlotlyJS.savefig(s1,"results/scatter_prob/scatter_$var-$vec_string_plot-s1.pdf");
#     PlotlyJS.savefig(s1,"results/scatter_prob/scatter_$var-$vec_string_plot-s1.eps");
# end