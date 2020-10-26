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
mu = 0.0000000000001;  # Resource mean
zeta = 2; # Resource variability scaling
edensity = 4000; # Resource energy density kJ/gram
mass = 10; # KILOGRAMS
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 50; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 10.0;
configurations = 100000;

cgut, cfat = acrossdaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,kmax,tmax_bout,cyears,configurations)

R"plot($cfat,type='l')"

mean(cfat)

zetavec = collect(1.0:0.001:2.0)
timedeath = Array{Int64}(undef,length(zetavec));
meanfat = Array{Float64}(undef,length(zetavec));
cvfat = Array{Float64}(undef,length(zetavec));

for i=1:length(zetavec)
    zeta = zetavec[i];
    alpha = 2;
    cgut, cfat = acrossdaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,gut_type,kmax,tmax_bout,cyears,configurations)

    meanfat[i] = mean(cfat);
    cvfat[i] = std(cfat)/mean(cfat);
    
    td = findall(x->x<1,cfat);
    if length(td) == 0
        timedeath[i] = cyears*365;
    else
        timedeath[i] = td[1];
    end
end
