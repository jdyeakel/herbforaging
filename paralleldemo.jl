


using Distributed
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using LinearAlgebra


l = 100;
A = SharedArray{Float64}(10,10);

ivec = repeat(collect(1:10),outer=10);
jvec = repeat(collect(1:10),inner=10);
parvec = [ivec jvec];

@time @sync @distributed for k=1:l
    i = parvec[k,1];
    j = parvec[k,2];
    A[i,j] = rand() + 2*rand() - 3*rand() + rand()^2;
end
