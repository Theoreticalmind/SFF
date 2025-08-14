using Plots
using DelimitedFiles
using Printf

function analyze_and_plot_sff()
    L=20
    #V = 0.0
    
    filename ="/Users/photon/Downloads/my_projects_jl/sff/sff_data_L20_Np4_V0.0_D0.1_NR10_b0.5_s0.01_200.dat"
    println("Analyzing data from: $filename")

    if !isfile(filename)
        println("Error: Data file not found!")
        return
    end

    data = readdlm(filename, ',', Float64, skipstart=1)
    time_data = data[:, 1]
    sff_data = data[:, 2]
    println("Successfully loaded $(length(time_data)) data points.")

    # Analyze the Plateau
    plateau_start_index = floor(Int, 0.90 * length(time_data))
    plateau_indices = plateau_start_index:length(time_data)

    plateau_value = mean(sff_data[plateau_indices])
    @printf("Plateau Position ≈ %.6e\n", plateau_value)

    #  Analyze the Ramp Slope    
    t_ramp_min = 10.0
    t_ramp_max = 200.0
    ramp_indices = findall(t -> t_ramp_min <= t <= t_ramp_max, time_data)

    if isempty(ramp_indices)
        println("Warning: Ramp region not found in the time window [$(t_ramp_min), $(t_ramp_max)]")
    else
        # Select the data for linear fitting
        t_ramp = time_data[ramp_indices]
        sff_ramp = sff_data[ramp_indices]


        log_t = log10.(t_ramp)
        log_sff = log10.(sff_ramp)

        X = hcat(ones(length(log_t)), log_t)
        
        coeffs = X \ log_sff
        intercept = coeffs[1]
        slope = coeffs[2]

        @printf("Ramp Slope ≈ %.4f\n", slope)

        # Generate the fitted line for plotting
        sff_fit = 10 .^ (intercept .+ slope .* log_t)
    end

    p = plot(time_data, sff_data,
        label="SFF Data", lw=2.5, color=:red,
        xaxis=:log10, yaxis=:log10,
        xlabel="Time (t)", ylabel="SFF F(t)",
        #title="SFF Analysis (L=$L, Np=$Np, V=$V)",
        legend=:topleft, framestyle=:box,
        size=(500, 400))

    hline!(p, [plateau_value],
        label="Plateau ≈ $(@sprintf("%.2e", plateau_value))",
        ls=:dashdot, lw=2, color=:green)

    if @isdefined(sff_fit)
        plot!(p, t_ramp,  sff_fit,
            label="Ramp Fit (Slope ≈ $(@sprintf("%.3f", slope)))",
            ls=:dot, lw=2.5, color=:blue)
    end

    println("\nDisplaying annotated plot.")
    #savefig("sff_plot_ED_L$(L)_V$(V).pdf")
    display(p)
end

analyze_and_plot_sff()