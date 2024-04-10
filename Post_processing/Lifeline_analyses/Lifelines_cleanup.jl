## Code to uniformize the lifelines in order to process them later on

# should make sure that each file is equally long (1000 columns) (start at 3000 s) and end at 4000 script


# load packages
using DelimitedFiles
using Random
# cd("U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run1_DynField/Lifelines")

# Lifeline_number = 42116
writing_switch_lifeline_files = 0                                                     # 1 to overwrite the data
writing_switch_sorted_files = 0
number_of_lifelines = 4000

function Lifeline_cleaner(Lifeline_number::Int64)
    cd("U:/MicroSynC/Lars/PhD/Paper_scale_down/FLUENT/DPM/Run1-DynField_5-gl/Lifelines")
    try Lifeline_data = readdlm("lifeline_$Lifeline_number.out",',');

        # make sure that it starts at 3000.1 s
        starting_time = 3000.1;                                             # moment at which DPM method started
        dt = 0.1;                                                           # timestep in DPM method
        end_time = 4000;

        if Lifeline_data[1,1] > starting_time                               # if particle was injected later on
            LL_injection_time = Lifeline_data[1,1];                         # retrieve data point of injection
            Added_steps = Int(round((LL_injection_time - starting_time)/dt));           # how many rows should be added with NaN data
            N_Nan_columns = size(Lifeline_data,2)-1;                                                     

            tvec = Array(starting_time:dt:(LL_injection_time-dt));          # time vector to be added
            NaN_Matrix = zeros(Added_steps, N_Nan_columns)*NaN;
            Front_added_matrix =  [tvec NaN_Matrix];
            Lifeline_data = [Front_added_matrix ; Lifeline_data];
        end  
        
        if Lifeline_data[end,1] < end_time                             # to recover the mistake... now we remove all data smaller than 4000 at the end of the lifeline
            timediff = diff(Lifeline_data[:,1])
            remove_rows_start = argmin(timediff)                                   # from begin to here should be kept
            
            Lifeline_data = Lifeline_data[1:remove_rows_start, :];
        end

        if writing_switch_lifeline_files == 1
            # open("Lifeline_$Lifeline_number.out","w") do io
            open("Lifeline_datawrite_try.out","w") do io
            writedlm(io, Lifeline_data,',') end
            # rm("Lifeline_datawrite_try.out")
        end

        if writing_switch_sorted_files == 1

            for i in 1:length(vars)
                var = vars[i];
                var_data = [Lifeline_number ; Lifeline_data[:,i+1]];

                outfile = "sorted_lifelines/movie-$var.out"
                    open(outfile, "a") do io
                writedlm(io,var_data')
                end
            end
            # open("Lifeline_$Lifeline_number.out","w") do io
            # open("sorted_lifelines/CO.out","w") do io
            # writedlm(io, clCO',',') end
            # rm("Lifeline_datawrite_try.out") 
        end
    catch
        println("Lifeline_$Lifeline_number.out does not exist")
    end    
    # return var_data
end


#       For the analysis script, use the histcount from matlab (https://stackoverflow.com/questions/70606493/matlab-histcounts-in-julia)

starting_number = 0;
# lifeline 158888 did not exist somehow
lifeline_array = Array(starting_number:1:number_of_lifelines::Int);
vars = ["x", "y", "z", "clCO","clH2","clCO2","qCO","qH2","qCO2"]

@time begin
for Lifeline_number in lifeline_array
    Lifeline_cleaner(Lifeline_number)
    if mod(Lifeline_number, 1000) == 0
        println("Lifeline $Lifeline_number written")
    end
   # print("Wrote lifeline $Lifeline_number in file")
end
println("All lifelines written")

end

