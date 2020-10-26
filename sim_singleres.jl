using Distributed
using SharedArrays
using Distributions
using RCall
using LinearAlgebra


include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/trait_and_rate_functions.jl");
include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/withindaysim_singleres.jl");
include("/home/jdyeakel/Dropbox/PostDoc/2020_herbforaging/src/acrossdaysim_singleres.jl");

# TESTRUN
rho = 1;
alpha = 2; # Resource dispersion
mu = 0.00000000001;  # Resource mean
zeta = 1; # Resource variability scaling
edensity = 4000; # Resource energy density kJ/gram
mass = 10; # KILOGRAMS
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 1000; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 0.5;
configurations = 100000;

gr, cgut, cfat, ginfo,gprob = acrossdaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,kmax,tmax_bout,cyears,configurations);

R"plot($cfat,type='l')"

R"plot($ginfo,$gprob,type='b')"

mean(cfat)
# Integrated energetic state / total possible
rfit = sum(cfat)/(maximum(cfat)*cyears*365)



# SIMULATE ACROSS ZETA
reps = 10;
zetavec = collect(1.0:0.001:2.0);
rfit = Array{Float64}(undef,length(zetavec));
rho = 1;
alpha = 2; # Resource dispersion
mu = 0.00000000001;  # Resource mean
zeta = 1; # Resource variability scaling
edensity = 4000; # Resource energy density kJ/gram
mass = 10; # KILOGRAMS
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 1000; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 1;
configurations = 100000;

for i=1:length(zetavec)
    zeta = zetavec[i];

    rfitvec = Array{Float64}(undef,reps)
    for r = 1:reps
        gr, cgut, cfat, ginfo,gprob = acrossdaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,kmax,tmax_bout,cyears,configurations);
        rfitvec[r] = sum(cfat)/(maximum(cfat)*cyears*365)
    end

    # meanfat[i] = mean(cfat);
    # cvfat[i] = std(cfat)/mean(cfat);
    rfit[i] = mean(rfitvec);

    percentdone = floor((i/length(zetavec))*100);
    if mod(percentdone,10) == 0
        println(percentdone)
    end
end

R"plot($rfit,pch=16)"