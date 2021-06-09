using Distributed
using Distributions
using LinearAlgebra
using RCall
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/trait_and_rate_functions.jl");
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/smartpath.jl");
@everywhere include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/withindaysim_singleres.jl");

#NOTE: run this over mass with alpha = 2.1

alphavec = collect(2.1:0.1:4);
zetavec = collect(1:0.01:2);

mcalc = Array{Float64}(undef,length(alphavec),length(zetavec));
varcalc = Array{Float64}(undef,length(alphavec),length(zetavec));
manalytical = Array{Float64}(undef,length(alphavec),length(zetavec));
varanalytical = Array{Float64}(undef,length(alphavec),length(zetavec));

for i = 1:length(alphavec)
    for j = 1:length(zetavec)
        #Resource
        rho = 10.0;
        mu = 0.00000000001;  # Resource mean
        alpha = alphavec[i];
        zeta = zetavec[j];
        #Consumer (kg)
        mass = 10^5;

        # Draw distance to next resource
        # First draw encounter rate
        draws = 100000;

        #Consumer population density: individuals/m^2
        n = indperarea(mass);
        #bite size
        beta = bite_size_allo(mass); # mass in kg, bite size in g/bite

        # mu*(1/beta) = resource density = bites/m^2
        # rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
        m = rho*mu*(1/beta);

        #Adjusted for competitive landscape
        mprime = m/n;
        alphaprime = alpha*n^(zeta-2);

        #Define Gamma Distribution for resource availability
        gammadist = Gamma(alphaprime,mprime/alphaprime); #mean = alpha * m/alpha

        rg = rand(gammadist,draws);
        distance_to_resource = rand.(Exponential.(1.0 ./ rg));

        mcalc[i,j] = mean(distance_to_resource);
        varcalc[i,j] = var(distance_to_resource);

        #Analytical
        manalytical[i,j] = (alpha*n^(zeta+1))/(m*(alpha*n^zeta-n^2));
        varanalytical[i,j] = (alpha^3*n^(3*zeta-4))/(m^2*(alpha*n^(zeta-2)-2)*(alpha*n^(zeta-2)-1)^2);

    end
end

cvanalytical = sqrt.(varanalytical)./manalytical;


namespace = smartpath("figures/sim_v_pred_variance.pdf")
R"""
pdf($namespace,width=6,height=5)
plot($(varcalc),$(varanalytical),log='xy',xlab='Simulated distance variance ',ylab='Predicted distance variance')
dev.off()
"""

namespace = smartpath("figures/distancesd_alpha_zeta.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$alphavec,y=$zetavec,z=sqrt($(varcalc))/1000)
dev.off()
"""

namespace = smartpath("figures/distancemean_alpha_zeta.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$alphavec,y=$zetavec,z=($(mcalc)/1000))
dev.off()
"""


namespace = smartpath("figures/distanceCV_alpha_zeta.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$alphavec,y=$zetavec,z=($(cvanalytical)))
dev.off()
"""



#Assuming alpha = 100 (i.e. large), calculate mean(distance) and var(distance) as a function of MASS and ZETA

zetavec = collect(1:0.1:2);
massexpvec = collect(0:0.1:5);
massvec = 10 .^massexpvec; #kg
rho = 5.0*10^-9;
alpha = 3; # Resource dispersion
mu = 1;  # Resource mean
foragehours = 2;
tmax_bout = foragehours*60*60; # Set at 1/2 day (6) hours

DMarray = Array{Float64}(undef,length(zetavec),length(massvec));
DVarray = Array{Float64}(undef,length(zetavec),length(massvec));

DMsimarray = Array{Float64}(undef,length(zetavec),length(massvec));
DVsimarray = Array{Float64}(undef,length(zetavec),length(massvec));

for i=1:length(zetavec)
    for j=1:length(massvec)

        zeta = zetavec[i];
        mass = massvec[j];


        velocity = find_velocity(mass); # m/s
        beta = bite_size_allo(mass); # mass in kg, bite size in g/bite
        #Consumer population density: individuals/m^2
        ndensity = indperarea(mass); #individuals/m^2
        forageseconds = copy(tmax_bout); #seconds
        homerangediameter = velocity*forageseconds; #meters
        homerangearea = pi*(homerangediameter/2)^2; #meters^2
        # n = ndensity*homerangearea; #inds/homerange
        n = ndensity;
        # n=1
        # rho * mu * (1/beta) = bite encounters/m^2 ~ bite encounter rate
        m = rho*mu*(1/beta);
        #Adjusted for competitive landscape
        mprime = m/n;
        alphaprime = alpha*(n^(zeta-2));

        #analytical
        mdist = alphaprime/(mprime*(alphaprime - 1));
        vardist = alphaprime^3/((mprime^2)*((alphaprime - 1)^2)*(alphaprime -2));

        DMarray[i,j] = mdist;
        DVarray[i,j] = vardist;

        draws = 10000;
        gammadist = Gamma(alphaprime,mprime/alphaprime); #mean = alpha * m/alpha
        rg = rand(gammadist,draws);
        distance_to_resource = rand.(Exponential.(1.0 ./ rg));


        # gprob, ginfo, tout = withindaysim_singleres(rho,alpha,mu,zeta,edensity,mass,teeth,kmax,tmax_bout,configurations);
        #This is a hack:
        # gprob[findall(x->isnan(x)==true,gprob)].=0;
        DMsimarray[i,j] = mean(distance_to_resource);
        DVsimarray[i,j] = var(distance_to_resource);
    end
end

namespace = smartpath("figures/distancemean_mass_zeta_ndensity.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$zetavec,y=$massexpvec,z=log($(DMarray)),xlab='zeta',ylab='Mass 10^i (kg)')
dev.off()
"""

namespace = smartpath("figures/distancevar_mass_zeta_ndensity.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$zetavec,y=$massexpvec,z=log($(DVarray)),xlab='zeta',ylab='Mass 10^i (kg)')
dev.off()
"""

namespace = smartpath("figures/distancemeanSIM_mass_zeta.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$zetavec,y=$massexpvec,z=log($(DMsimarray)),xlab='zeta',ylab='Mass 10^i (kg)',zlim=c(1,10))
dev.off()
"""

namespace = smartpath("figures/distancevarSIM_mass_zeta.pdf")
R"""
library(fields)
pdf($namespace,width=7,height=7)
image.plot(x=$zetavec,y=$massexpvec,z=log($(DVsimarray)),xlab='zeta',ylab='Mass 10^i (kg)',zlim=c(1,10))
dev.off()
"""
