function withindaysim_singleres(
    rho,
    alpha, # Resource dispersion
    mu,  # Resource mean
    zeta, # Resource variability scaling
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    teeth, # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
    kmax, # 50 in Sevilleta NOTE: I *think* controls bin size?
    tmax_bout,
    configurations
    )

    probability = Array{Float64}(undef,kmax+1);
    kinfo = Array{Float64}(undef,kmax+1);
    # probability = SharedArray{Float64}(kmax+1);
    # kinfo = SharedArray{Float64}(kmax+1);
    data = zeros(Float64,configurations);


    let tictoc = 0, bad_draw = 1
        #Rebuild if there is a bad draw
        while bad_draw == 1

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
            ndensity = indperarea(mass); #individuals/m^2
            # n=1; #Need to think harder about how n -> mprime

            forageseconds = copy(tmax_bout); #seconds
            homerangediameter = velocity*forageseconds; #meters
            homerangearea = pi*(homerangediameter/2)^2; #meters^2
            n = ndensity*homerangearea; #inds/area


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
            bites = 0;
            
            probability = probability .* 0.0;
            kinfo = kinfo .* 0.0;
            data = data .* 0.0;

            # Slow down organism if they have more choices
            # modvelocity = maximum([tweight[target],1/num_res])*velocity;

            # bitesperday = Array{Int64}(undef,configurations);

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
                # BITES = Array{Int64}(undef,1);
                let number_of_successes = 0, nearest_resource = 0, t=0.0, distance_to_resource = 0.0, gut = 0.0, tic = 0
                
                    while t < tmax_bout
                        tic += 1


                        # Draw distance to next resource
                        # First draw encounter rate
                        rg = rand(gammadist);
                        distance_to_resource = rand(Exponential(1.0/rg));
                            
                        #The forager will move towards the closest resource
                        ttravel = distance_to_resource/velocity; # distance / velcity = time
                        
                        # to avoid zero food encounters, FORCE the first or second draw a shorter distance!
                        if tic == 1 || tic == 2
                            if ttravel > (tmax_bout*0.9);
                                ttravel = tmax_bout/10;
                            end
                        end
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
                        if (t + ttravel) < tmax_bout
                            number_of_successes += 1;
                            
                            # Time passes while chewing
                            t += tchew; #time
                            t_chew += tchew; #time

                            # Pass Mouth-Unit to Gut Total (boundary conditions in across-day sim)
                            # resgain is kJ per volume (mouth volume)
                            gut += beta; #grams/bite
                            bites += 1;
                            
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
                    # BITES[1] = bites;
                end
                        
                # total_kilojoules=dot((resgain),number_of_successes); #.*epsilon
                # avg_digestibility = epsilon .* (number_of_successes/sum(number_of_successes));
                # bitesperday[config] = BITES[1];
                data[config] = GUT[1] * edensity; #grams * kJ/gram = kJ returns
            end
            # mbitesperday = mean(bitesperday);
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

            # if any(isnan.(probability)) == true
            #     tictoc += 1;
            #     bad_draw = 1;
            #     if tictoc == 10
            #         println("Problem building within day return distribution")
            #         break
            #     end
            # else
            #     bad_draw = 0;
            # end
            bad_draw = 0;
        end #end while
    end #end let

    # Calculate average times
    t_travel /= configurations;
    t_chew /= configurations;
    # t_wait /= configurations;
    bites /= configurations;
    
    # Can tuples be > 2?
    tout = tuple(t_travel,t_chew,bites);
    
    
    #note: probability matrix should have dims (kmax)
    return (probability[1:kmax-1], kinfo[1:kmax-1], tout);
end
