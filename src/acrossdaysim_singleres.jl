function acrossdaysim_singleres(
    rho,
    alpha, # Resource dispersion
    mu,  # Resource mean
    zeta, # Resource variability scaling
    edensity, # Resource energy density kJ/gram
    mass, # KILOGRAMS
    teeth, # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
    gut_type,
    kmax, # 50 in Sevilleta NOTE: I *think* controls bin size?
    tmax_bout, # Set at 12 hours (43200 seconds)
    cyears, #consumer years
    configurations)

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
    greturnprob, greturninfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
    
    # Visually check distributions
    greturnmean = dot(greturnprob,greturninfo);
    # R"barplot($greturnprob)";

    
    # CALCULATE twait
    # Maximum gut capacity
    maxgut = gut_volume_g(mass, gut_type) * edensity; # grams * kJ/gram = kJ

    # NOTE THIS IS TOO LONG
    mrt = mean_retention_time(mass, gut_type); # seconds / PARTICLE

    particle_mass = mean_particle_mass(mass, gut_type); #gram / particle

    # Passage rate of food (rate of flow from gut to body)
    passrate = (1/mrt) * particle_mass; #particle/s * gram / particle = grams/s

    # Single day gut passage
    gutpass = passrate * secday * edensity; # grams/s * s/day * kJ/gram = kJ/day

    # Metabolic loss per day
    maxfatstorage, cost_basalrate, cost_fieldrate = find_metabolism(mass);

    #Active time
    activehrs = 12;
    field_cost = cost_fieldrate * (activehrs*60*60);    #kJ/s * seconds = KJ
    rest_cost = cost_basalrate * ((24 - activehrs)*60*60); #KJ/s * seconds = KJ

    # How long does it take to empty gut from maxgut to 1/2 maxgut?
    # twait = (maxgut/2)*(1/passrate); #grams * s/grams = seconds

     #### NEEDS WORK ####
    # Digestibility
    # NOTE What is the relationship between retention time and digestibility????
    # Closer to 1 for ruminants
    # Less than 1 for simple digestive systems

    # Because MRT is positively linked to the digestive efficiency of a herbivore (Foose, 1982; Udén and Van Soest, 1982; Clauss et al., 2007b),
    # epsilon = (max_retention_time(mass) - mean_retention_time(mass, gut_type))/(max_retention_time(mass) - (min_retention_time(mass)));
    epsilon = 0.5;

    daysincyears = Int64(floor(365*cyears));  

    cgut = zeros(daysincyears);
    cfat = zeros(daysincyears);

    #Starting conditions
    cgut[1] = maxgut;
    cfat[1] = maxfatstorage;

    probline = cumsum(greturnprob);
    
    for t=2:daysincyears
        if cfat[t-1] > 0.0
            #Draw daily return
            fdraw = rand();
            gutreturn = greturninfo[findall(x->x>fdraw,probline)[1]];

            #Problem... gut passage rate has to be multiplied by gut contents!

            deltagut = (gutreturn - cgut[t-1]*gutpass);
            cgut[t] = minimum([maximum([cgut[t-1] + deltagut,0.0]),maxgut]);

            deltafat = (epsilon*cgut[t-1]*gutpass - field_cost - rest_cost);
            cfat[t] = minimum([maximum([cfat[t-1] + deltafat,0.0]),maxfatstorage]);
        end
    end

    return cgut, cfat


end