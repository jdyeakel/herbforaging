function richness_mass_eval(reps,rhovec,massexpvec,zvec,anat_traits,res_traits)

    teeth = anat_traits[:teeth];
    gut_type = anat_traits[:gut_type];

    mu = res_traits[:mu];
    alpha = res_traits[:alpha];
    edensity = res_traits[:edensity];


    # zvec = [1,2];
    # rhovec = collect(0.010:0.01:0.500); #0.04; #used with ndensity*homerange
    # rhovec = collect(1:2:50).*10^-9; #used with ndensity
    # massexpvec = collect(0:0.2:4.4); #exponent for mass in KG
    # reps = 10;
    pv1 = [repeat(collect(1:length(massexpvec)),inner = length(rhovec)) repeat(collect(1:length(rhovec)),outer=length(massexpvec))];
    pv2 = [repeat(collect(1:reps),inner=length(rhovec)*length(massexpvec)) repeat(pv1,outer=reps)];
    its = size(pv2)[1];

    mpropalive = Array{Float64}(undef,length(zvec),length(massexpvec),length(rhovec));
    mcovfatres = Array{Float64}(undef,length(zvec),length(massexpvec),length(rhovec));
    mcovrelfatres = Array{Float64}(undef,length(zvec),length(massexpvec),length(rhovec));
    for j = 1:length(zvec)
        zeta = zvec[j];
        # alpha = 5; # Resource dispersion
        # mu = 1;  # Resource mean
        # edensity = 18.2; # Resource energy density kJ/gram (from YKR)
        # teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
        # gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
        kmax = 100; # 50 in Sevilleta NOTE: I *think* controls bin size?
        foragehours = 2;
        tmax_bout = foragehours*60*60; # Set at 1/2 day (6) hours
        configurations = 10000;

        propalive = SharedArray{Float64}(reps,length(massexpvec),length(rhovec));
        covfatres = SharedArray{Float64}(reps,length(massexpvec),length(rhovec));
        covrelfatres = SharedArray{Float64}(reps,length(massexpvec),length(rhovec)); #NOTE: CV of relative values not good

        @time @sync @distributed for i = 1:its
            #parameter position
            p = pv2[i,:];
            #parameter value
            r = p[1];
            mass = 10^massexpvec[p[2]];
            rho = rhovec[p[3]];
            
            cyears = 5.08*mass^0.35; #Mean expected lifespan (yrs); Calder 1984; 10;

            gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
            #This is a hack:
            gprob[findall(x->isnan(x)==true,gprob)].=0;
            # R"plot($ginfo,$gprob,type='b')"
            gr, cgut, cfat, ctraits = acrossdaysim_singleres(gprob,ginfo,edensity,mass,gut_type,cyears);
            relfatres = cfat/ctraits[1];
            days_dead = findall(iszero,relfatres);
            if length(days_dead) > 0
                days_alive = Float64(days_dead[1]);
            else
                days_alive = (cyears * 365);
            end

            propalive[p[1],p[2],p[3]] =  days_alive / (cyears * 365);
            covfatres[p[1],p[2],p[3]] = sqrt(var(cfat[1:Int64(floor(days_alive))]))/mean(cfat[1:Int64(floor(days_alive))]);
            covrelfatres[p[1],p[2],p[3]] = sqrt(var(relfatres[1:Int64(floor(days_alive))]))/mean(relfatres[1:Int64(floor(days_alive))]);
        end
            
        mpropalive[j,:,:] = mean(propalive,dims=1)[1,:,:];
        mcovfatres[j,:,:] = mean(covfatres,dims=1)[1,:,:];
        mcovrelfatres[j,:,:] = mean(covrelfatres,dims=1)[1,:,:];

    end

    return mpropalive, mcovfatres, mcovrelfatres

end