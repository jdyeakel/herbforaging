if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/herbforaging/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/loadfuncs.jl");
end




# using Distributed
# @everywhere using SharedArrays
# @everywhere using Distributions
# @everywhere using RCall
# @everywhere using LinearAlgebra

# @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/trait_and_rate_functions.jl");
# @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/withindaysim_singleres.jl");
# @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/acrossdaysim_singleres.jl");
# @everywhere include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/smartpath.jl");


# SIMULATE ACROSS ZETA
reps = 50;
zetavec = collect(1.0:0.01:2.0);
massexpvec = collect(1:0.1:6);
rfit = SharedArray{Float64}(length(zetavec)*length(massexpvec));
rho = 0.1;
alpha = 2; # Resource dispersion
mu = 0.00000000001;  # Resource mean
zeta = 1; # Resource variability scaling
edensity = 4000; # Resource energy density kJ/gram
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 1000; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 1;
configurations = 100000;

paramvalues = [repeat(zetavec,outer=length(massexpvec)) repeat(massexpvec,inner=length(zetavec))];


@time @sync @distributed for i=1:(size(paramvalues)[1])
    
    zeta = copy(paramvalues[i,1]);
    massexp = copy(paramvalues[i,2]);
    mass = 10^massexp;

    # Within-parameter repetitions
    rfitvec = Array{Float64}(undef,reps);

    gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
    #This is a hack:
    gprob[findall(x->isnan(x)==true,gprob)].=0;

    for r = 1:reps
        gr, cgut, cfat = acrossdaysim_singleres(gprob,ginfo,edensity,mass,gut_type,cyears);
        gprob[findall(x->isnan(x)==true,gprob)].=0;
        rfitvec[r] = sum(cfat)/(maximum(cfat)*cyears*365)
    end

    # meanfat[i] = mean(cfat);
    # cvfat[i] = std(cfat)/mean(cfat);
    rfit[i] = mean(rfitvec);

    # percentdone = floor((i/length(zetavec))*100);
    # if mod(percentdone,10) == 0
    #     println(percentdone)
    # end
end

rfitmatrix = reshape(rfit,(length(zetavec),length(massexpvec)));

namespace = smartpath("figures/fitness_zetamass.pdf")
R"""
pdf($namespace,width=6,height=5)
image(x=$zetavec,y=$massexpvec,z=log($rfitmatrix),pch=16,xlab='Zeta',ylab='Mass')
dev.off()
"""

