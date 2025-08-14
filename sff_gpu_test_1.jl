using LinearAlgebra
using Random
using Distributions
using Plots
using Printf
using Combinatorics
using DelimitedFiles
using CUDA
using CUDA.CUSPARSE # For sparse matrices on the GPU

# This function now prepares the data for the Hamiltonian on the CPU
# and then assembles a sparse matrix on the GPU.
function build_interacting_hamiltonian_gpu(L::Int, Np::Int, t1::Float32, t2::Float32, V::Float32,
                                           wij_disorder::Dict{Tuple{Int,Int},Float32})

    basis = [Tuple(c) for c in combinations(1:L, Np)]
    basis_size = length(basis)
    state_to_index = Dict(state => i for (i, state) in enumerate(basis))

    # Using Coordinate (COO) format for sparse matrix construction
    I = Int32[]
    J = Int32[]
    V_vals = Float32[]

    for i in 1:basis_size
        state = basis[i]

        # --- Diagonal: Interaction term ---
        interaction_energy = 0.0f0
        for site_j in 1:(L-1)
            if (site_j in state) && ((site_j + 1) in state)
                interaction_energy -= V
            end
        end
        if interaction_energy != 0
            push!(I, i)
            push!(J, i)
            push!(V_vals, interaction_energy)
        end

        # --- Off-diagonal: Hopping terms (t1, t2) ---
        for site_k in state
            neighbor = site_k + 1
            if neighbor <= L && !(neighbor in state)
                hopping_t = (site_k % 2 == 1) ? t1 : t2
                new_state_tuple = Tuple(sort!([s for s in state if s != site_k] ∪ [neighbor]))
                j = state_to_index[new_state_tuple]

                push!(I, i)
                push!(J, j)
                push!(V_vals, -hopping_t)
            end
        end

        # --- Off-diagonal: Disorder hopping (w_ij) ---
        for (sites, w_ij) in wij_disorder
            site_a, site_b = sites
            if isodd(site_a) == isodd(site_b) continue end

            if (site_a in state) && !(site_b in state)
                new_state_tuple = Tuple(sort!([s for s in state if s != site_a] ∪ [site_b]))
                j = state_to_index[new_state_tuple]
                push!(I, i)
                push!(J, j)
                push!(V_vals, -w_ij)
            end
        end
    end

    # Create a sparse matrix directly on the GPU
    H_sparse_gpu = CuSparseMatrixCSR(cu(I), cu(J), cu(V_vals), (basis_size, basis_size))
    # Symmetrize the matrix since we only added one triangle
    return H_sparse_gpu + H_sparse_gpu'
end

#= This function is now adapted to work with CuArrays
function calculate_sff_gpu(eigenvalues::CuVector{Float32}, time_array::CuVector{Float32}, β::Float32)::Vector{Float32}
    E = reshape(eigenvalues, :, 1)
    t = reshape(time_array, 1, :)

    # Perform calculations on the GPU
    Z_plus = sum(exp.(-(β .+ 1im .* t) .* E), dims=1)
    Z_minus = sum(exp.(-(β .- 1im .* t) .* E), dims=1)
    Z_beta = sum(exp.(-β .* E))

    if Z_beta ≈ 0
        return zeros(Float32, length(time_array))
    end

    sff_val_gpu = vec((Z_plus .* Z_minus) ./ (Z_beta^2))

    # Move the final result back to the CPU
    return real(Array(sff_val_gpu))
end=#



# This function is now adapted to work with CuArrays and is numerically stable.
function calculate_sff_gpu(eigenvalues::CuVector{Float32}, time_array::CuVector{Float32}, β::Float32)::Vector{Float32}
    # --- FIX: Stabilize the calculation by shifting the energy spectrum ---
    # Find the ground state energy (the minimum eigenvalue).
    E_min = minimum(eigenvalues)
    # Shift all energies so that the ground state is at 0. This is physically
    # equivalent and prevents numerical overflow in the exp() function.
    E_shifted = eigenvalues .- E_min

    # Reshape for broadcasting, using the shifted energies
    E = reshape(E_shifted, :, 1)
    t = reshape(time_array, 1, :)

    # Perform calculations on the GPU
    Z_plus = sum(exp.(-(β .+ 1im .* t) .* E), dims=1)
    Z_minus = sum(exp.(-(β .- 1im .* t) .* E), dims=1)
    Z_beta = sum(exp.(-β .* E))

    # This check is still good practice, in case all eigenvalues were identical
    # and the temperature was infinite (beta=0), though unlikely here.
    if Z_beta ≈ 0
        return zeros(Float32, length(time_array))
    end

    sff_val_gpu = vec((Z_plus .* Z_minus) ./ (Z_beta^2))

    # Move the final result back to the CPU
    return real(Array(sff_val_gpu))
end

function main()
    # System parameters (L=18, Np=6 as requested)
    L = 18
    Np = 6
    V = 2.0f0      # Use Float32 for GPU performance
    τ = 2.0f0
    Δ = -0.1f0
    t1_top = (τ +  Δ) / 2.0f0
    t2_top = (τ -  Δ) / 2.0f0

    NR = 10
    σ = 0.01f0

    N_samples = 100
    time_points = 500
    # Move time_array to the GPU once
    time_array_d = cu(exp10.(range(-1, 4, length=time_points)))
    β = 0.5f0

    # Calculate basis size for info messages
    basis_size = binomial(L, Np)
    println("--- CUDA-based Exact Diagonalization SFF Simulation ---")
    println("System Size L=$L, Particle Number Np=$Np, Interaction V=$V")
    println("Disorder Strength σ=$σ in NR=$NR sites")
    println("Many-body Hilbert space dimension: $basis_size x $basis_size")
    println("Averaging over $N_samples disorder samples...")
    println("Using device: $(CUDA.name(CuDevice(0)))")

    # Array to store results on the CPU
    sff_top_samples = zeros(Float32, length(time_array_d), N_samples)

    r_start = div(L - NR, 2) + 1
    r_end = r_start + NR - 1
    R_sites = r_start:r_end
    
    # How many eigenvalues to compute. For SFF, you typically need all of them.
    # NOTE: Calculating all eigenvalues for a 18564x18564 matrix is very demanding.
    # We'll try to get them all, but this may fail on GPUs with less memory.
    # If it fails, reduce `num_eigenvalues`.
    num_eigenvalues = basis_size - 2 # svds requires k < min(n, m)

    for i in 1:N_samples
        if i % 10 == 0; @printf("  Sample %d/%d\n", i, N_samples); end

        # Generate disorder on the CPU
        wij_disorder = Dict{Tuple{Int,Int},Float32}()
        for site_i in R_sites, site_j in R_sites
            if site_i < site_j && isodd(site_i + site_j)
                wij_disorder[(site_i, site_j)] = rand(Normal(0, σ / sqrt(NR)))
            end
        end

        # Build Hamiltonian on the GPU as a sparse matrix
        H_top_gpu = build_interacting_hamiltonian_gpu(L, Np, t1_top, t2_top, V, wij_disorder)

        # Calculate eigenvalues on the GPU
        # Using svds for sparse matrices. 'sr' means smallest real part.
        # Note: svds computes singular values, which are sqrt(eigvals) for H'H.
        # For a Hermitian matrix, eigvals(H'H) = eigvals(H)^2.
        # This is an approximation. For true eigenvalues of a sparse matrix,
        # you might need a library like KrylovKit.jl with GPU support.
        # For simplicity, we form the dense matrix here if it fits in memory.
        try
            H_dense_gpu = CuMatrix(H_top_gpu)
            evals_top_gpu = eigvals(Hermitian(H_dense_gpu))
            sff_top_samples[:, i] = calculate_sff_gpu(sort(evals_top_gpu), time_array_d, β)
        catch e
            if isa(e, CUDA.OutOfMemoryError)
                println("\nCUDA out of memory when creating the dense Hamiltonian.")
                println("Consider using a smaller L or a GPU with more memory.")
                return # Exit if memory runs out
            else
                rethrow(e)
            end
        end

        # It's good practice to free GPU memory explicitly when done with large arrays in a loop
        CUDA.reclaim()
    end

    avg_sff_top = vec(mean(sff_top_samples, dims=2))

    println("Saving data to file...")
    filename = "sff_data_L$(L)_Np$(Np)_V$(V)_D$(Δ)_NR$(NR)_b$(β)_s$(σ)_GPU.dat"
    data_to_save = hcat(Array(time_array_d), avg_sff_top)

    open(filename, "w") do io
        writedlm(io, ["Time" "SFF"], ',')
        writedlm(io, data_to_save, ',')
    end
    println("Data successfully saved to $filename")

end

# Check if a CUDA-capable GPU is available before running
if CUDA.functional()
    main()
else
    println("No functional CUDA device found. Please ensure you have a")
    println("CUDA-capable GPU and the necessary drivers installed.")
end