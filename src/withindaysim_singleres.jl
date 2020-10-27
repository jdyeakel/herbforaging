function withindaysim_singleres(
    rho,
    alpha, # Resource dispersion
    mu,  # Resource mean
    zeta, # Resource variability scaling
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    teeth, # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
    kmax, # 50 in Sevilleta NOTE: I *think* controls bin size?
    tmax_bout, # Set at 12 hours (43200 seconds)
    configurations
    )


    ######################
    # CONSUMERS
    ######################

    #Grams per bite
    # from shipley 94 "the scaling of intake rate"
    # mass in kg
    # bite size in g
    # Why elephants have trunks Pretorius 2015
    # bite_grams = (0.002 * (mass)^0.969); # grams
    # bite_rate = 0.37 * mass^(-0.024); #(bites/s)

    velocity = find_velocity(mass); # m/s

    # CALCULATE tchew
    # NOTE: BITE SIZE SEEMS SMALL
    # bite_rate = bite_rate_allo(mass); # mass in kg, bite/s
    beta = bite_size_allo(mass); # mass in kg, bite size in g/bite

    # mouthrate = bite_rate * bite_size; # bite/s * g/bite = grams/s
    # 1/mouthrate is seconds/1 gram

    chewrate = chew_allo(mass, teeth); #g/s
    tchewgram = 1/chewrate; #s/g
    tchew = tchewgram * beta; #s/g * g/bite = s/bite

    
   


    ###################
    # RESOURCES
    ###################

    # SCALE ENCOUNTERS TO BITE! (i.e. each encounter should be a 'bite')
    # NOTE: BUILD IN ZETA SCALING!

    #Consumer population density: individuals/m^2
    n = indperarea(mass);

    # mu*(1/beta) = resource density = bites/m^2
    # rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
    m = rho*mu*(1/beta);

    #Adjusted for competitive landscape
    mprime = m/n;
    alphaprime = alpha*n^(zeta-2);

    #Define Gamma Distribution for resource availability
    gammadist = Gamma(alphaprime,mprime/alphaprime); #mean = alpha * m/alpha
    
    t_travel = 0.0;
    t_chew = 0.0;
    
    # probability = Array{Float64}(undef,kmax+1).*0.0;
    # kinfo = Array{Float64}(undef,kmax+1).*0.0;
    probability = SharedArray{Float64}(kmax+1);
    kinfo = SharedArray{Float64}(kmax+1);
    data = zeros(Float64,configurations);

    # Slow down organism if they have more choices
    # modvelocity = maximum([tweight[target],1/num_res])*velocity;
    
    # @sync @distributed 
    for config = 1:configurations
        # Counters :: within configuration
        # number_of_successes = 0; # Each success is a bite
        # nearest_resource = 0;
        # t=0.0;
        # distance_to_resource = 0.0;
        
        # # Energetic Returns!
        # gut = 0.0;
        GUT = Array{Float64}(undef,1);
        let number_of_successes = 0, nearest_resource = 0, t=0.0, distance_to_resource = 0.0, gut = 0.0
        
            while t < tmax_bout

                # Draw distance to next resource
                distance_to_resource = rand(Exponential(1.0/rand(gammadist)));
                    
                #The forager will move towards the closest resource
                ttravel = distance_to_resource/velocity; # distance / velcity = time
                t += ttravel; # time
                t_travel += ttravel; # time

                # # Digest while traveling
                # # What has been digested?
                # digested = gut*passrate*ttravel; #grams * 1/time * time = grams
                # # Subtract this from the gut
                # gut -= digested; # grams
                # gut = maximum([gut,0.0]); # grams
                # # Add this to the ereturns modified by digestibility epsilon
                # ereturns += epsilon*digested; # grams
                
                # If the forager successfully reaches the resource
                # Each resource is equal to 1 mouthful
                if tmax_bout > (distance_to_resource/velocity)
                    number_of_successes += 1;
                    
                    # Time passes while chewing
                    t += tchew; #time
                    t_chew += tchew; #time

                    # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
                    # resgain is kJ per volume (mouth volume)
                    gut += beta; #grams/bite
                    
                    
                    # gut = minimum([gut,maxgut]); #grams

                    # # digest while chewing
                    # # What has been digested?
                    # digested = gut*passrate*tchew; # grams * 1/time * time = grams
                    # # Subtract this from the gut
                    # gut -= digested; #grams
                    # gut = maximum([gut,0.0]); #grams
                    # # Add this to the ereturns modified by digestibility epsilon
                    # ereturns += epsilon*digested; #grams
                end
                
                # # If you fill your stomach, wait until it is at xcapacity before starting again
                # if gut == maxgut
                #     t += twait; #time 
                #     t_wait += twait; #time
                #     digested = gut*passrate*twait; #grams * 1/time * time = grams
                #     gut -= digested; #grams
                #     gut = maximum([gut,0.0]); #grams
                #     ereturns += epsilon*digested; #grams
                # end 
            end
            GUT[1] = gut;
        end
                
        # total_kilojoules=dot((resgain),number_of_successes); #.*epsilon
        # avg_digestibility = epsilon .* (number_of_successes/sum(number_of_successes));
        
        data[config] = GUT[1] * edensity; #grams * kJ/gram = kJ returns
    end
    
    datamin = minimum(data);
    datamax = maximum(data);
    if datamin!=datamax
        for config = 1:configurations
            probability[convert(Int64,trunc(kmax*(data[config]-datamin)/(datamax-datamin)))+1]+=1.0/configurations;
        end
        k = collect(datamin:(datamax-datamin)/kmax:datamax);
        for i = 1:kmax
            kinfo[i]=k[i];
        end
    else
        probability[1]=1.;
        k=datamin;
        kinfo[1]=datamin;
    end

    # Calculate average times
    t_travel /= configurations;
    t_chew /= configurations;
    # t_wait /= configurations;
    
    # Can tuples be > 2?
    tout = tuple(t_travel,t_chew);
    
    
    #note: probability matrix should have dims (kmax)
    return (probability, kinfo, tout);
end
