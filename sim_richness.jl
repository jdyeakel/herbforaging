if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/herbforaging/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2020_herbforaging/src/loadfuncs.jl");
end


zvec = [1,2];
rhovec = collect(1:5:100).*10^-9; #used with ndensity
massexpvec = collect(0:0.2:4.4);
reps = 10;
res_traits = (mu = 1, alpha = 3, edensity = 18.2);

teethvec = ["bunodont","acute/obtuse lophs", "lophs and non-flat", "lophs and flat"];
gut_typevec = ["caecum", "colon", "non-rumen foregut", "rumen foregut"];
its = length(teethvec)*length(gut_typevec);

mpa = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(zvec),length(massexpvec),length(rhovec));
mcovf = Array{Float64}(undef,length(teethvec),length(gut_typevec),length(zvec),length(massexpvec),length(rhovec));
let ijtic = 0
    for i=1:length(teethvec)
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

filename = "data/richness/herbivoretraits.jld"
namespace = smartpath(filename)
@save namespace reps rhovec massexpvec zvec teethvec gut_typevec res_traits its mpa mcovf


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
        image($(massexpvec),$rhovec,$((mcovf[i,j,1,:,:])),main=$(string(teethvec[i],"; ",gut_typevec[j])),xlab='Body mass 10^i',ylab = 'Richness')
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
        image($(massexpvec),$rhovec,$((mcovf[i,j,2,:,:])),main=$(string(teethvec[i],"; ",gut_typevec[j])),xlab='Body mass 10^i',ylab = 'Richness')
        """
    end
end
R"dev.off()"

#Find min richness for body mass 
minrho = Array{Float64}(undef,length(massexpvec));
for i=1:
for j=1:length(massexpvec)
    impossiblerhos = findall(x->x < 1., mpa[1,1,1,j,:]);
    if length(impossiblerhos) > 0
        minrho[j] = rhovec[last(impossiblerhos)];
    else
        minrho[j] = maximum(rhovec);
    end
end
lineplot(massexpvec,minrho)