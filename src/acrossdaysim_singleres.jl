function acrossdaysim_singleres(
    gprob, 
    ginfo,
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    gut_type,
    cyears) #consumer years

    # rho = 1;
    # alpha = 2; # Resource dispersion
    # mu = 0.00000001;  # Resource mean
    # zeta = 1; # Resource variability scaling
    # edensity = 4000; # Resource energy density kJ/gram
    # mass = 10; # KILOGRAMS
    # teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
    # kmax = 50; # 50 in Sevilleta NOTE: I *think* controls bin size?
    # tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
    # configurations = 100000;

    #Second in a day
    secday = 60*60*24;

    # Daily Gut Return Distributions
    #first build within-day gut returns distributions (kJ)
    # gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);

    # If there is a bad within-day simulation, rerun until there is not!
    # let tictoc = 0
    #     while any(isnan.(gprob)) == true
    #         tictoc += 1;
    #         gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
    #         if tictoc == 10
    #             println("Problem within day")
    #         end
    #     end
    # end
    
    # Visually check distributions
    greturnmean = dot(gprob,ginfo);
    # R"plot($ginfo,$gprob,type='b')"

    
    # CALCULATE twait
    # Maximum gut capacity
    maxgut = gut_volume_g(mass, gut_type) * edensity; # grams * kJ/gram = kJ

    # NOTE THIS IS TOO LONG
    mrt = mean_retention_time(mass, gut_type); # seconds / PARTICLE

    particle_mass = mean_particle_mass(mass, gut_type); #gram / particle

    # Passage rate of food (rate of flow from gut to body)
    passrate = (1/mrt) * particle_mass; #particle/s * gram / particle = grams/s

    # Single day gut passage
    # passage rate of a kJ within a single particle (needs to be multiplied by kJ in stomach to get total kJ passing in a day)
    gutpass = passrate * secday * edensity; # grams/s * s/day * kJ/gram = kJ/day

    # Metabolic loss per day
    activehrs = 12;
    maxfatstorage, cost_basalrate, cost_fieldrate = find_metabolism(mass,activehrs);
    #maxfatstorage = kJ
    #costs = kJ/day

    #Active time
    field_cost = cost_fieldrate;    #kJ/day
    rest_cost = cost_basalrate; #KJ/day

    # How long does it take to empty gut from maxgut to 1/2 maxgut?
    # twait = (maxgut/2)*(1/passrate); #grams * s/grams = seconds

     #### NEEDS WORK ####
    # Digestibility
    # NOTE What is the relationship between retention time and digestibility????
    # Closer to 1 for ruminants
    # Less than 1 for simple digestive systems

    # Because MRT is positively linked to the digestive efficiency of a herbivore (Foose, 1982; Udén and Van Soest, 1982; Clauss et al., 2007b),
    # epsilon = (max_retention_time(mass) - mean_retention_time(mass, gut_type))/(max_retention_time(mass) - (min_retention_time(mass)));
    epsilon = 0.1;

    daysincyears = Int64(floor(365*cyears));  

    cgut = zeros(daysincyears);
    cfat = zeros(daysincyears);
    gr = zeros(daysincyears);

    #Starting conditions
    cgut[1] = maxgut;
    cfat[1] = maxfatstorage;

    probline = cumsum(gprob);
    probline = probline/maximum(probline);
    
    for t=2:daysincyears


        if cfat[t-1] > 0.0

             # We could build in environmental stress here
            # (Prevalence of Good <-> Bad days w/ autocorrelation)
            gooddaydraw = rand();

            if gooddaydraw < 1
                #good day
                #Draw daily return
                fdraw = rand();
                #draw from greturns (allowed to be > stomach size)
                gutreturndraw = findall(x->x>fdraw,probline); #kJ
                gutreturn = ginfo[gutreturndraw[1]];
            else
                #bad day
                gutreturn = 0.0;
            end

            # if length(gutreturndraw) == 0
                # fdraw = rand();
                #draw from greturns (allowed to be > stomach size)
                # gutreturndraw = ginfo[findall(x->x>fdraw,probline)[1]]; 
                #kJ
                #if fdraw is very large, nothing on probline will be > fdraw
                # gr[t] = maximum(ginfo);
            # else
            gr[t] = gutreturn;
            # end

            #Change in stomach contents
            deltagut = (gr[t] - cgut[t-1]*gutpass);
            cgut[t] = minimum([maximum([cgut[t-1] + deltagut,0.0]),maxgut]);

            deltafat = (epsilon*cgut[t-1]*gutpass - field_cost - rest_cost);

            #We might want an incorporation compartment... where gut contents are passed to a compartment, and added to fat at a constant rate...

            cfat[t] = minimum([maximum([cfat[t-1] + deltafat,0.0]),maxfatstorage]);
        end
    end
        
    ctraits = tuple(maxfatstorage)

    return gr, cgut, cfat, ctraits


end