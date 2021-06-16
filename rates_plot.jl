if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/herbforaging/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/loadfuncs.jl");
end

massexpvec = collect(0:0.2:4.4);

teethvec = ["bunodont","acute/obtuse lophs", "lophs and non-flat", "lophs and flat"];
gut_typevec = ["caecum", "colon", "non-rumen foregut", "rumen foregut"];
res_traits = (mu = 1, alpha = 3, edensity = 18.2);


filename = "figures/fig_retentiontime.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=5,height=4)
plot($massexpvec,$(mean_retention_time.(10 .^massexpvec,gut_typevec[1])./(60*60*24)),type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='Gut retention time (days)',log='y')
"""
for i=2:length(gut_typevec)
    R"""
    lines($massexpvec,$(mean_retention_time.(10 .^massexpvec,gut_typevec[i])./(60*60*24)),col=pal[$i],lwd=2)
    """
end
R"""
legend(0,10,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


passrate = Array{Float64}(undef,length(gut_typevec),length(massexpvec));
for i=1:length(gut_typevec);
    #NOTE: Put these together:
    maxgut = gut_volume_g.(10 .^massexpvec, gut_typevec[i]) .* res_traits[:edensity]; # grams * kJ/gram = kJ
    # NOTE THIS IS TOO LONG
    mrt = mean_retention_time.(10 .^massexpvec, gut_typevec[i]); # seconds / PARTICLE
    particle_mass = mean_particle_mass.(10 .^massexpvec, gut_typevec[i]); #gram / particle
    # Passage rate of food (rate of flow from gut to body)
    passrate[i,:] = (1 ./ mrt) .* particle_mass .* ((60*60*24)); #particle/s * gram / particle = grams/s
end

filename = "figures/fig_passagerates.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=5,height=4)
plot($massexpvec,$(passrate[1,:]),type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='Passage rate (grams/day)',log='y',ylim=c(min($passrate),max($passrate)))
"""
for i=2:length(gut_typevec)
    R"""
    lines($massexpvec,$(passrate[i,:]),col=pal[$i],lwd=2)
    """
end
R"""
legend(0,10000,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


# NOTE: PLOT BITES PER AREA FOR FULL RANGE OF RHOVEC