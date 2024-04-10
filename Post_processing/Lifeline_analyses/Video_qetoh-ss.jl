using JLD2
using StatsBase
using LinearAlgebra
using PlotlyJS
using LaTeXStrings
using WebIO
using FileIO
using DelimitedFiles

function CO_uptake(cCO)
    # inputs are concentrations in mM
    qcomax = 1.459;                 # molco/cmol/h
    K_m = 0.042;                    # mM
    K_i = 0.25;                     # mM^2
        
    q = qcomax*cCO/(K_m+cCO+(cCO^2)/K_i); # mol/molx/h
    return q
end

function H2_uptake(cH2, cCO)
    # inputs are concentrations in mM
    qh2max = 2.565;                 # molH2/cmol/h
    K_m = 0.025;                    # mol/m3
    K_i = 0.025;                    # mol/m3
    
    q = qh2max*cH2/(K_m+cH2)*( 1 /(1 + (cCO/K_i))); # mol/molx/h
    return q
end


function fun_qEtOH(clCO, clH2)
    qH2 = H2_uptake(clH2, clCO)
    qCO = CO_uptake(clCO)
    qEtOH = 1/6*qCO + 1/6*qH2

    return qEtOH
end

function load_coords(PATH, num_lifelines)
    # f = jldopen("$PATH/coords_movie_$num_lifelines.jld2");
    # xcoords = f["x_data"]'
    # ycoords = f["y_data"]'
    # zcoords = f["z_data"]'
    # clCO = f["clCO_data"]'
    # clH2 = f["clH2_data"]'
    
    f_x_coords = readdlm("$PATH/movie-x.out");
    f_y_coords = readdlm("$PATH/movie-y.out");
    f_z_coords = readdlm("$PATH/movie-z.out");
    f_clCO = readdlm("$PATH/movie-clCO.out");
    f_clH2 = readdlm("$PATH/movie-clH2.out");

    xcoords = f_x_coords[ [1:601 ; 603:1689 ; 1691:num_lifelines+3 ],2:end];
    ycoords = f_y_coords[  [1:601 ; 603:1689 ; 1691:num_lifelines+3 ],2:end];
    zcoords = f_z_coords[ [1:601 ;603:1689 ; 1691:num_lifelines+3 ],2:end];
    clCO = f_clCO[ [1:601 ; 603:1689 ; 1691:num_lifelines+3 ],2:end];
    clH2 = f_clH2[ [1:601 ; 603:1689 ; 1691:num_lifelines+3 ],2:end];

    # xcoords[602,:] = [];

    # calculate EtOH uptake rate
    qEtOH = zeros(size(clCO))

    # calculate qEtOH
    for i = 1:size(clCO,1)
        for j = 1:size(clCO,2)
            qEtOH[i,j] = fun_qEtOH(clCO[i,j], clH2[i,j])
            if (qEtOH[i,j ] < 0)
                qEtOH[i,j] = 0;
            elseif qEtOH[i,j] > 1
                qEtOH[i,j] = 1;
            end
        end
    end

    xcoords[1,:] .= -5;
    ycoords[1,:] .= -5;
    zcoords[1,:] .= 0;
    qEtOH[1,:] .= 1e-9;
    clCO[1,:] .= 1e-6;

    return xcoords, ycoords, zcoords, clCO, qEtOH
end

function plotfun(xcoords, ycoords, zcoords, clCO, qEtOH_norm)
    trace = cat(scatter(x = xcoords, y = ycoords, z = zcoords,
    mode="markers",
    marker=attr(
        size = 1 .+ qEtOH_norm.*24,
        color=clCO, showscale=true, line_width=0, colorbar_title_text = "<i>   c<sub>L,CO</sub> </i> <br> (mol m<sup>-3</sup>)",
        colorbar = attr(len = 0.6, tickfont_size = 14, title_font_size = 16),
        colorscale="Viridis", colorbar_exponentformat="power", cmin=0, cmax = 0.25,      # set color to an array/list of desired values
        opacity=1), type="scatter3d"), dims = 1);

    return trace
end


function fLayout(time)
     # src = base64encode(image)
    # img = load("$PATH/Animation/reactor-darker.png");
    t = round(time;  digits=1);

    layout = Layout(
        title=attr(x = 0.7, y = 0.88, text = "<i>t</i> = $t s"), 
        scene_xaxis=attr(title = "", range=[-5, 5], showgrid=false, gridwidth=0.3, gridcolor="grey", 
        zeroline=false,zerolinewidth=1, zerolinecolor="grey", type="linear", showline=false,
        ticks="", tickwidth=2, tickcolor="grey", ticklen=5, position = 0, anchor = 0, automargin=true, visible = 0, tickvals = []),
        scene_yaxis=attr(title = "", range=[-5,10],showgrid=false, gridwidth=0.3, gridcolor="grey",
        zeroline=false, showline=false, zerolinewidth=1, zerolinecolor="grey", 
        type="linear", ticks="", tickwidth=2, tickcolor="grey", ticklen=5, 
        automargin=true, anchor = 0, position = 0, visible = 0, tickvals = []),
        scene_zaxis=attr(title = "", range=[0,24.5],showgrid=false, gridwidth=0.3, gridcolor="grey",
        zeroline=false, showline=false, zerolinewidth=1, zerolinecolor="grey", 
        type="linear", ticks="", tickwidth=2, tickcolor="grey", ticklen=5, 
        automargin=true, anchor = 0, position = 0, visible = 0, tickvals = []),
        # legend_title="Legend Title",
        font=attr(
        family="times",
        size=12,
        color="black"),
        paper_bgcolor="rgb(255,255,255)",
        plot_bgcolor="rgb(255,255,255)",
        width = 600, height = 900, showlegend=false,
        # updatemenus = updatemenus,
        # sliders = sliders,
        aspectmode="manual", aspectratio=attr(x=1, y=2, z=6),
        scene_camera = attr(
            up=attr(x=0, y=0, z=1),                 # z-axis points upwards
            center=attr(x=0, y=0, z=0),
            eye=attr(x=2, y=2, z=2),
            projection = "orthographic"), 
        scene_bgcolor = "rgb(255,255,255)",
        scene_xaxis_backgroundcolor = "rgba(0,0,0,0)",
        scene_yaxis_backgroundcolor = "rgba(0,0,0,0)",
        scene_zaxis_backgroundcolor = "rgba(0,0,0,0)",
        # <a href="https://freeimage.host/i/ydFmYP"><img src="https://iili.io/ydFmYP.md.png" alt="ydFmYP.md.png" border="0"></a>        
        # <a href="https://freeimage.host/i/Hox77Tb"><img src="https://iili.io/Hox77Tb.md.png" alt="Hox77Tb.md.png" border="0"></a>
        # <a href="https://freeimage.host/i/Hoxchue"><img src="https://iili.io/Hoxchue.md.jpg" alt="Hoxchue.md.jpg" border="0"></a>
      #  <a href="https://freeimage.host/i/Hoxl3P4"><img src="https://iili.io/Hoxl3P4.md.png" alt="Hoxl3P4.md.png" border="0"></a>
        images = [attr(source="https://iili.io/Hoxl3P4.md.png", 
        xref="paper",
        yref="paper",
        x=-0.025,
        y=1.02,
        sizex=1.0,
        sizey=0.9,
        sizing="fill",
        opacity=0.5)] 
        )
    return layout
end

function fPlotHist(qEtOH_norm)
    vector = qEtOH_norm[:];
    bins = 0:0.1:1;
    h1 = fit(Histogram, vector, bins)
    h = normalize(h1, mode=:probability)

    max_axis_value = ceil(maximum(h.weights), sigdigits=:1)
    layout = Layout(
        title="",
        xaxis_title="<i>q<sub>EtOH</sub>/q<sub>EtOH,max</sub> </i>",
        xaxis=attr( showgrid=true, gridwidth=0.3, gridcolor="lightgrey", 
            zeroline=true,zerolinewidth=1, zerolinecolor="black", type="linear",
            ticks="inside", tickwidth=2, tickcolor="black", ticklen=5),
        yaxis_title="probability",
        yaxis=attr(range=[0, max_axis_value], showgrid=true, gridwidth=0.3, gridcolor="lightgrey",
            zeroline=true, showline=true, zerolinewidth=1, zerolinecolor="black", 
            ticks="inside", tickwidth=2, tickcolor="black", ticklen=5, 
            automargin=true, anchor = 0, position = 0),
        legend_title="Legend Title",
        font=attr(
            family="times",
            size=16,
            color="black"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width = 700, height = 500 );

    trace = bar(x = Array(h.edges[1]), y = h.weights, marker=attr(color = "royalblue"));
    hplot = plot(trace, layout)
    return hplot
end


# tstart, dt, tend = 3000.2, 0.1, 3650;
num_lifelines = 2000;
# start_time_0 = 100;                     # from the first frame 
# start_time = 305.4;                     # for the current batch of frames to be made (can also start at start_time_0)
# end_time = 400;                         # for the current batch of frames to be made
# plot_freq = 0.1;                        # every second (or 0.1 s when = 0.1)
# tvec = tstart:dt:tend;    

movie_maker = 1; save_fig_test = 1; histogram_qEtOH = 0;
dt = 0.1;
start_time_index = 1999;                     # from the first frame 
end_time_index = 4002;                     # from the first frame 

t = 0:dt:(1000-dt);
start_time_0 = t[start_time_index];                     # from the first frame 
start_time = start_time_0;                     # for the current batch of frames to be made (can also start at start_time_0)
end_time = t[end_time_index];                         # for the current batch of frames to be made
plot_freq = dt;       
tvec = start_time:dt:end_time;              


PATH = "U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run1-DynField_5-gl/Lifelines/sorted_lifelines";

#       Get data from data_pkvals script
@time begin
xcoords, ycoords, zcoords, clCO, qEtOH  = load_coords(PATH, num_lifelines)
end
# xcoords, ycoords, zcoords, clCO, qEtOH  = xcoords', ycoords', zcoords', clCO', qEtOH';
println("data loaded!")

# normalize it based upon its maximum uptake rate
# qEtOH_norm = minimum(qEtOH, 1)
qEtOH_norm = qEtOH/1;


time_index_st = round(Int, start_time_0/dt)+1;
x_pos_start, y_pos_start, z_pos_start = xcoords[:,time_index_st], ycoords[:, time_index_st], zcoords[:, time_index_st];
clCO_start_time, qEtOH_norm_start_time = clCO[:,time_index_st], qEtOH_norm[:,time_index_st];

if histogram_qEtOH == 1
    h = fPlotHist(qEtOH_norm);
    display(h)
end


trace = plotfun(x_pos_start, y_pos_start, z_pos_start, clCO_start_time, qEtOH_norm_start_time);
pl = plot(trace, fLayout(start_time))
frame_index = round(Int,(start_time)/plot_freq+1 );
println("frame-index = $frame_index")

if save_fig_test == 1
    savefig(pl, "$PATH/frames2/frame_$frame_index.png",  width=600, height=900, scale=5)
    # savefig(pl, "$PATH/Animation_qEtOH/frame0.pdf")
end
display(pl)

# make animation
n_frames = floor(Int,length(start_time:dt:end_time)/(plot_freq/dt) + 1);      # number of frames to be made
frames  = Vector{PlotlyFrame}(undef, n_frames)
global frame_iter = 1;

if movie_maker == 1
    updatemenus = [attr(type="buttons", 
        active=0, y=1, x=1.1, #(x,y) button position 
        buttons=[attr(label="Play", method="animate", args=[nothing, attr(frame=attr(duration=5, redraw=true), 
        transition=attr(duration=0), fromcurrent=true, mode="immediate")])])];

    sliders = [attr(active=1, minorticklen=0, steps=[attr(label="$((k-1)*plot_freq+start_time)", method="animate",
                args=[["$k"], # match the frame[:name]
                    attr(mode="immediate", transition=attr(duration=0),
                        frame=attr(duration=5, redraw=true))
                    ]) for k in 1:n_frames ] )];    

    println("begin looping over time")

    # loop over the whole time in steps of 1 s
    for time = start_time:dt:end_time
        time_index = round(Int,time/dt) + 1;
        
        x_pos, y_pos, z_pos = xcoords[:,time_index], ycoords[:, time_index], zcoords[:, time_index];
        clCO_time, qEtOH_norm_time = clCO[:,time_index], qEtOH_norm[:,time_index];
            
        # # save as frame
        frame_index = round(Int,(time)/plot_freq+1 );
        # frames[frame_index] = frame(data=[scatter(x=x_pos, y=y_pos, z=z_pos,
        #     mode="markers",
        #     marker=attr(size=clCO_pks[:,time_index].*5, color=clCO[:,time_index]), 
        #     type="scatter3d")],
        #     layout=attr(title_text="Time = $time s"),            #update title
        #     name="$frame_index", #frame name; it is passed to slider 
        #     traces=[0]) # this means that the above data update the first trace (here the unique one) 

        # make a new plot 
        trace_new = plotfun(x_pos, y_pos, z_pos, clCO_time, qEtOH_norm_time)
        pl = plot(trace_new, fLayout(time))
        sleep(0.1)
        savefig(pl, "$PATH/frames2/frame$(frame_index).png",  width=600, height=900, scale=5)

        
        if mod(time_index, 10) == 0
            println("time = $time s")
        end
 
    end
end

# make a video using ffmpeg in the BASH 
#   1)      set home directory to be the $PATH/Animation folder
#   2)      set the ffmpeg command
#               ffmpeg -framerate 25 -start_number 1000 -i frames2/frame%d.png -vcodec libx264 -s 600x900 Video/frames_25fps.mp4
#                         frames per s    input         encoding       resol      output 

# then, to make it compatible with powerpoint use the command
#   3)      ffmpeg -i Video/frames_10fps.mp4 -c:v libx264 -preset slow -profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 15 Video/frames_10fps_ppt.mp4
#                           import              codex       ..              needed                              needed?   15 is high quality
#
