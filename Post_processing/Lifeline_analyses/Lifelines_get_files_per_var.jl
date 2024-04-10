
# makes datasets for every variable as a function of t: x, y, z, cCO, cH2, cCO2, qCO, qH2, qCO2
using CSV
using DelimitedFiles
using DataFrames
using StringEncodings

# LIST OF VARIABLES IN YOUR .OUT FILE PER LIFELINE
vars = ["x", "y", "z", "clCO","clH2","clCO2","qCO","qH2","qCO2"];
# PATH TO YOUR DIRECTORY WITH LIFELINE.OUT FILES
PATH = "U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run5-DynField_10-gl/Lifelines";
writing_switch = 1
max_lifelines = 41714;      # take from folder with lifelines


function Lifeline_sorter(Lifeline_number::Int64, PATH, writing_switch)
    cd("$PATH")
    try Lifeline_data = readdlm("Lifelines/lifeline_$Lifeline_number.out",',');

        t = Lifeline_data[:,1];
    
        # check whether upcoming row is correct, specify end point...   
        i = 1; ind_keep = []; ind_remove = []; tf_3650 = 0;
        while i <= length(t)-1
            if t[i] == 3650;  tf_3650 = 1; end
            if t[i+1] < t[i] && tf_3650 == 0;
                ind_remove  = push!(ind_remove, i+1);
                ind_keep = push!(ind_keep,i);
                i += 1;
            elseif tf_3650 == 1 && t[i] == 3650
                ind_keep = push!(ind_keep,i);
            elseif tf_3650 == 1
                ind_remove  = push!(ind_remove, i);
            else
                ind_keep = push!(ind_keep,i);
            end
            i += 1;
        end


        if writing_switch == 1

            for i in 1:length(vars)
                var = vars[i];
                var_data = Lifeline_data[ind_keep,i+1];

                outfile = "sorted_lifelines/$var.out"
                    open(outfile, "a") do io
                writedlm(io,var_data')
                end
            end
        end
        return t, Lifeline_data
    catch
        println("Lifeline_$Lifeline_number.out does not exist")
    end
end

@time begin
for Lifeline_number in Array(0:max_lifelines)
    Lifeline_sorter(Lifeline_number, PATH, writing_switch)
    if Lifeline_number%100 == 0
        println("Lifeline $Lifeline_number written")
    end
    end
end

println("Wrote all lifelines --- finished!!")

