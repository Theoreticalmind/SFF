using LinearAlgebra
using Combinatorics
using Printf

function build_hamiltonian(L::Int, Np::Int, t1::Float64, t2::Float64, V::Float64)
    # This function remains unchanged.
    basis = [Tuple(c) for c in combinations(1:L, Np)]
    basis_size = length(basis)
    
    state_to_index = Dict(state => i for (i, state) in enumerate(basis))

    H_many_body = zeros(Float64, basis_size, basis_size)

    for i in 1:basis_size
        state = basis[i] 

        interaction_energy = 0.0
        
        for site_j in 1:(L-1)
            if (site_j in state) && ((site_j + 1) in state)
                interaction_energy += V 
            end
        end
        H_many_body[i, i] = interaction_energy

        for site_k in state
            neighbor = site_k + 1 
            if neighbor <= L && !(neighbor in state) 

                hopping_t = isodd(site_k) ? t1 : t2
                
                new_state_tuple = Tuple(sort!([s for s in state if s != site_k] ∪ [neighbor]))
                
                j = state_to_index[new_state_tuple]
                
                H_many_body[j, i] = -hopping_t
                H_many_body[i, j] = -hopping_t
            end
        end
    end

    return Hermitian(H_many_body)
end


function main()
    L = 20         
    Np = 4          
    V = 2.0         
    
    τ = 2.0        
    Δ = 0.5        
    t1 = (τ + Δ) / 2.0
    t2 = (τ - Δ) / 2.0

    basis_size = binomial(L, Np)

    println("--- Exact Diagonalization Simulation ---")
    println("System Size L=$L, Particle Number Np=$Np")
    println("Hopping t1=$t1, t2=$t2, Interaction V=$V")
    println("Many-body Hilbert space dimension: $basis_size x $basis_size")
    
    H = build_hamiltonian(L, Np, t1, t2, V)
    
    println("\nCalculating eigenvalues...")
    eigenvalues = eigvals(H)
    println("Calculation complete.")

    println("\n--- Results ---")

    @printf "All %d Eigenvalues:\n" length(eigenvalues)
    for i in 1:length(eigenvalues)
        @printf "  E_%d = %.6f\n" (i-1) eigenvalues[i]
    end

    output_filename = "/home/antpc/Downloads/eigenvalues_sq_20_4.txt"
    try
        open(output_filename, "w") do f
            println(f, "# System Size L=$L, Particle Number Np=$Np")
            println(f, "# Hopping t1=$t1, t2=$t2, Interaction V=$V")
            
            for E in eigenvalues
                println(f, E)
            end
        end
        println("\nSuccessfully saved all eigenvalues to '$output_filename'")
    catch e
        println("\nAn error occurred while writing to file: $e")
    end
    # --- MODIFICATION END ---
end

main()