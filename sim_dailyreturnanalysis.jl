using Distributed
using UnicodePlots
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using RCall
@everywhere using LinearAlgebra

@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/trait_and_rate_functions.jl");
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/withindaysim_singleres.jl");
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/acrossdaysim_singleres.jl");
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/smartpath.jl");

# TESTRUN
rho = 0.01;
alpha = 2; # Resource dispersion
mu = 1;  # Resource mean
zeta = 1; # Resource variability scaling
edensity = 18.2; # Resource energy density kJ/gram (from YKR)
mass = 10000; # KILOGRAMS
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 100; # 50 in Sevilleta NOTE: I *think* controls bin size?
foragehours = 2;
tmax_bout = foragehours*60*60; # Set at 1/2 day (6) hours
cyears = 10;
configurations = 100000;

gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
#This is a hack:
gprob[findall(x->isnan(x)==true,gprob)].=0;
R"plot($ginfo,$gprob,type='b')"


gr, cgut, cfat = acrossdaysim_singleres(gprob,ginfo,edensity,mass,gut_type,cyears);
R"plot($cfat,type='l')"


# Possible fitness measure
# Integrated energetic state / total possible
rfit = sum(cfat)/(maximum(cfat)*cyears*365)



# SIMULATE ACROSS ZETA
zetavec = collect(1.0:0.01:2.0);
rho = 10;
alpha = 2; # Resource dispersion
mu = 0.00000000001;  # Resource mean
zeta = 1; # Resource variability scaling
edensity = 4000; # Resource energy density kJ/gram
mass = 1000; # KILOGRAMS
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 1000; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 1;
configurations = 100000;

m_returns = Array{Float64}(undef,length(zetavec));
var_returns = Array{Float64}(undef,length(zetavec));

@time @sync @distributed for i=1:length(zetavec)
    zeta = zetavec[i];

    rfitvec = Array{Float64}(undef,reps);

    gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
    #This is a hack:
    gprob[findall(x->isnan(x)==true,gprob)].=0;

    #Calculate the daily return expectation
    exp_return = dot(ginfo,gprob);
    exp_squarereturn = dot(ginfo .^2,gprob);
    m_returns[i] = exp_return;

    #Calculate the daily return variance
    var_returns[i] = exp_squarereturn - exp_return^2;

    #Build data (for testing)
    #NOTE: above variance equation works
    # global x = Array{Float64}(undef,0);
    # for k=1:length(ginfo)
    #     global x = [x;repeat([ginfo[k]],Int64(floor(gprob[k]*1000)))];
    # end

    # meanfat[i] = mean(cfat);
    # cvfat[i] = std(cfat)/mean(cfat);
    # rfit[i] = mean(rfitvec);

    percentdone = floor((i/length(zetavec))*100);
    if mod(percentdone,10) == 0
        println(percentdone)
    end
end

namespace = smartpath("figures/returns_mean_sd_zeta.pdf")
R"""
pdf($namespace,width=10,height=10)
par(mfrow=c(2,2))
plot($zetavec,$m_returns,pch=16,xlab='Zeta',ylab='Expected daily returns')
plot($zetavec,sqrt($var_returns),pch=16,xlab='Zeta',ylab='SD daily returns')
plot($zetavec,sqrt($var_returns)/$m_returns,pch=16,xlab='Zeta',ylab='CV daily returns')
dev.off()
"""

#Compare CV across mass, zeta

# SIMULATE ACROSS ZETA
zetavec = collect(1.0:0.01:2.0);
massexpvec = collect(0:0.1:5);
rho = 1;
alpha = 2.1; # Resource dispersion
mu = 0.00000000001;  # Resource mean
edensity = 4000; # Resource energy density kJ/gram
teeth = "bunodont"; # bunodont, acute/obtuse lophs, lophs and non-flat, lophs and flat
gut_type = "caecum"; # caecum, colon, non-rumen foregut, rumen foregut
kmax = 1000; # 50 in Sevilleta NOTE: I *think* controls bin size?
tmax_bout = 6*60*60; # Set at 1/2 day (6) hours (43200 seconds)
cyears = 1;
configurations = 100000;

paramvalues = [repeat(zetavec,outer=length(massexpvec)) repeat(massexpvec,inner=length(zetavec))];

m_returns = SharedArray{Float64}(length(zetavec)*length(massexpvec));
var_returns = SharedArray{Float64}(length(zetavec)*length(massexpvec));

@time @sync @distributed for i=1:(length(zetavec)*length(massexpvec))
    zeta = copy(paramvalues[i,1]);
    massexp = copy(paramvalues[i,2]);
    mass = 10^massexp;

    gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
    #This is a hack:
    gprob[findall(x->isnan(x)==true,gprob)].=0;

    #Calculate the daily return expectation
    exp_return = dot(ginfo,gprob);
    exp_squarereturn = dot(ginfo .^2,gprob);
    m_returns[i] = exp_return;

    #Calculate the daily return variance
    var_returns[i] = exp_squarereturn - exp_return^2;

    #Build data (for testing)
    #NOTE: above variance equation works
    # global x = Array{Float64}(undef,0);
    # for k=1:length(ginfo)
    #     global x = [x;repeat([ginfo[k]],Int64(floor(gprob[k]*1000)))];
    # end

    # meanfat[i] = mean(cfat);
    # cvfat[i] = std(cfat)/mean(cfat);
    # rfit[i] = mean(rfitvec);

    percentdone = (i/length(paramvalues))*100;
    if mod(percentdone,10) == 0
        println(percentdone)
    end
end

marray = reshape(m_returns,(length(zetavec),length(massexpvec)));
vararray = reshape(var_returns,(length(zetavec),length(massexpvec)));
cvarray = sqrt.(vararray) ./ marray;

namespace = smartpath("figures/cv_zeta_mass.pdf")
R"""
library(RColorBrewer)
pal=rev(colorRampPalette(brewer.pal(9,"Spectral"))($(length(massexpvec))));
pdf($namespace,width=7,height=5)
plot($zetavec,$(cvarray[:,1]),pch=16,col=pal[1],ylim=c(0,0.6))
"""
for i=2:length(massexpvec)
    R"points($zetavec,$(cvarray[:,i]),pch=16,col=pal[$i])"
end
R"""
legend(1,0.3,floor(10^$(massexpvec[[1,25,51]])),col=pal[c(1,25,51)],pch=16)
dev.off()
"""


