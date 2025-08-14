# SFF



in the main function you can see lots of parameter like:
    L = 20          
    Np = 4
    V = 2.0      
    τ = 2.0
    Δ = - 0.1
    t1_top = (τ +  Δ) / 2.0
    t2_top = (τ -  Δ) / 2.0
    


    #t1_triv = 1.0
    #t2_triv = 1.0  
    
    NR = 10        
    σ = 0.01        
    
    N_samples = 100   
    time_points = 500
    time_array = exp10.(range(-1, 4, length=time_points))
    β = 0.5


where L is no of sites, Np is no of fillings, V is interaction potential, \Delta is dimmerization paramter , \tau is hopping .


you only need to vary V and if possible on gpu just chnage the filling Np , chnage the \Delta negative and postive with the same value.
