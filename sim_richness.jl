if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/herbforaging/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/loadfuncs.jl");
end


zvec = [1,2];
# rhovec_orig = collect(0.1:1:100).*10^-9; #used with ndensity
rhoexpvec = collect(-11:0.03:-6)
rhovec = 10 .^ rhoexpvec

massexpvec = collect(0:0.1:4.4);
reps = 10;
res_traits = (mu = 1, alpha = 3, edensity = 18.2);

teethvec = ["bunodont","acute/obtuse lophs", "lophs and non-flat", "lophs and flat"];
gut_typevec = ["caecum", "colon", "non-rumen foregut", "rumen foregut"];
its = length(teethvec)*length(gut_typevec);

mpa = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(zvec),length(massexpvec),length(rhovec));
mcovf = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(zvec),length(massexpvec),length(rhovec));
let ijtic = 0
    @showprogress 1 "Computing..." for i=1:length(teethvec)
        for j=1:length(gut_typevec)

            anat_traits = (teeth=teethvec[i],gut_type=gut_typevec[j]);
            res_traits = (mu = 1, alpha = 3, edensity = 18.2);

            mpropalive, mcovfatres, mcovrelfatres = richness_mass_eval(reps,rhovec,massexpvec,zvec,anat_traits,res_traits);

            mpa[i,j,:,:,:] = mpropalive;
            mcovf[i,j,:,:,:] = mcovfatres;
            ijtic += 1;
            println(string("completed ",ijtic,"/",its))
        end
    end
end

filename = "data/richness/herbivoretraits2.jld"
namespace = smartpath(filename)
@save namespace reps rhovec massexpvec zvec teethvec gut_typevec res_traits its mpa mcovf
# @load namespace reps rhovec massexpvec zvec teethvec gut_typevec res_traits its mpa mcovf


filename = "figures/fig_richness_herbtraits_z1.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,width=12,height=12)
par(mfrow=c(4,4))
# image($(transpose(mcovf[1,1,1,:,:])))
"""
for i=1:length(teethvec)
    for j = 1:length(gut_typevec)
        R"""
        image($(massexpvec),$rhoexpvec,$((mcovf[i,j,1,:,:])),main=$(string(teethvec[i],"; ",gut_typevec[j])),xlab='Body mass 10^i',ylab = 'Richness')
        """
    end
end
R"dev.off()"


filename = "figures/fig_richness_herbtraits_z2.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,width=12,height=12)
par(mfrow=c(4,4))
"""
for i=1:length(teethvec)
    for j = 1:length(gut_typevec)
        R"""
        image($(massexpvec),$rhoexpvec,$((mcovf[i,j,1,:,:])),main=$(string(teethvec[i],"; ",gut_typevec[j])),xlab='Body mass 10^i',ylab = 'Richness')
        """
    end
end
R"dev.off()"



#Find min richness for body mass 
minrho = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(massexpvec));
for i=1:length(teethvec)
for j=1:length(gut_typevec)
for k=1:length(massexpvec)
    impossiblerhos = findall(x->x < 1., mpa[i,j,1,k,:]);
    if length(impossiblerhos) > 0
        minrho[i,j,k] = rhovec[last(impossiblerhos)];
    else
        minrho[i,j,k] = minimum(rhovec);
    end
end
end
end
minrhoz1 = minrho;
filename = "figures/fig_minrho_z1.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=6,height=5)
plot($massexpvec,$(minrho[1,1,:]),ylim=c(0,12*10^-8),type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='Minimal viable richness')
"""
for i=1:length(teethvec)
    for j=1:length(gut_typevec)
    R"lines($massexpvec,$(minrho[i,j,:]),col=pal[$j],lwd=2)"
    end
end
R"""
legend(3,12*10^-8,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


#Find min richness for body mass 
minrho = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(massexpvec));
for i=1:length(teethvec)
for j=1:length(gut_typevec)
for k=1:length(massexpvec)
    impossiblerhos = findall(x->x < 1., mpa[i,j,2,k,:]);
    if length(impossiblerhos) > 0
        minrho[i,j,k] = rhovec[last(impossiblerhos)];
    else
        minrho[i,j,k] = minimum(rhovec);
    end
end
end
end
minrhoz2 = minrho
filename = "figures/fig_minrho_z2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=6,height=5)
plot($massexpvec,$(minrho[1,1,:]),ylim=c(0,20*10^-8),type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='Minimal viable richness')
"""
for i=1:length(teethvec)
    for j=1:length(gut_typevec)
    R"lines($massexpvec,$(minrho[i,j,:]),col=pal[$j],lwd=2)"
    end
end
R"""
legend(3,20*10^-8,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


filename = "figures/fig_minrho_diff.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=6,height=5)
plot($massexpvec,$(minrhoz1[1,1,:]) - $(minrhoz2[1,1,:]),ylim=c(-8*10^-8,4*10^-8),type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='Diff z1-z2: Minimal viable richness')
"""
for i=1:length(teethvec)
    for j=1:length(gut_typevec)
    R"lines($massexpvec,$(minrhoz1[i,j,:]) - $(minrhoz2[i,j,:]),col=pal[$j],lwd=2)"
    end
end
R"""
lines($massexpvec,numeric(length($massexpvec)),lty=3)
legend(3,4*10^-8,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


pty = lineplot(massexpvec,minrho[1,:])
[lineplot!(pty,massexpvec,minrho[i,:]) for i=2:length(teethvec)]





# Derivatives!

filename = "figures/fig_dminrho_dmass_z1.pdf";
namespace = smartpath(filename);
drhomin_dmass = (-1) .* diff(minrhoz1[1,1,:]) ./ diff((10 .^ massexpvec));
R"""
library(RColorBrewer)
pal=brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=6,height=5)
plot($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmass,type='l',col=pal[1],lwd=2,xlab='Initial mass 10^i',ylab='(-) drhomin/dmass',xlim=c(0,1),ylim=c(0,12*10^-8))
"""
for i=1:length(teethvec)
    for j=1:length(gut_typevec)
        drhomin_dmass = (-1) .* diff(minrhoz1[i,j,:]) ./ diff((10 .^ massexpvec));
        R"lines($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmass,col=pal[$j],lwd=2)"
    end
end
# for i=1:length(teethvec)
#     for j=1:length(gut_typevec)
#         drhomin_dmass = (-1) .* diff(minrhoz2[i,j,:]) ./ diff((10 .^ massexpvec));
#         R"lines($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmass,col=pal[$j],lwd=2,lty=2)"
#     end
# end
R"""
legend(0.7,12*10^-8,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""

filename = "figures/fig_dminrho_dmass_diffz1z2.pdf";
namespace = smartpath(filename);
drhomin_dmassz1 = (-1) .* diff(minrhoz1[1,1,:]) ./ diff((10 .^ massexpvec));
drhomin_dmassz2 = (-1) .* diff(minrhoz2[1,1,:]) ./ diff((10 .^ massexpvec));
R"""
library(RColorBrewer)
pal=brewer.pal(length($gut_typevec),'Set1')
pdf($namespace,width=6,height=5)
plot($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmassz1-$drhomin_dmassz2,type='l',col=pal[1],lwd=2,xlab='Mass 10^i',ylab='(-) drhomin/dmass',xlim=c(0,1),ylim=c(-2*10^-8,2*10^-7.5))
"""
for i=1:length(teethvec)
    for j=1:length(gut_typevec)
        drhomin_dmassz1 = (-1) .* diff(minrhoz1[i,j,:]) ./ diff((10 .^ massexpvec));
        drhomin_dmassz2 = (-1) .* diff(minrhoz2[i,j,:]) ./ diff((10 .^ massexpvec));
        R"lines($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmassz1-$drhomin_dmassz2,col=pal[$j],lwd=2)"
    end
end
# for i=1:length(teethvec)
#     for j=1:length(gut_typevec)
#         drhomin_dmass = (-1) .* diff(minrhoz2[i,j,:]) ./ diff((10 .^ massexpvec));
#         R"lines($(massexpvec[1:(length(massexpvec)-1)]),$drhomin_dmass,col=pal[$j],lwd=2,lty=2)"
#     end
# end
R"""
lines(seq(0,5),numeric(6),lty=3)
legend(0.6,6*10^-8,$gut_typevec,col=pal,pch=16,cex=0.8)
dev.off()
"""


#ylim=c(0,2*10^-7.5),