using DrWatson
@quickactivate "S2DTMatrices"

using DelimitedFiles
using ARMLive
using Statistics

τ = parse(Int64, ARGS[4])
prefix = "$(ARGS[1])-$(ARGS[2])-$(ARGS[3])"
fnc = "$prefix-τ$τ-states.nc"
ds  = NCDataset(datadir("states","$(ARGS[5])",fnc))
dt  = ds["time"][:] .- Hour(6)
states = ds["states"][:]; nstates = maximum(states)
uniques = ds["unique"][:]
sizes = ds["sizes"][:] ./ 10
close(ds)

mkpath(datadir("transitionmatrix"))

statemat = zeros(nstates,nstates)

for it = 1 : (length(states)-1)
    x = states[it]
    y = states[it+1]
    if !iszero(x) && !iszero(y); statemat[x,y] += 1 end
end

statemat2 = statemat ./ sum(statemat,dims=2)

open(datadir("transitionmatrix","$prefix-τ$τ-transitionmatrix-basic.txt"),"w") do io
    write(io,"$(statemat2)")
end

for ii = 1 : nstates
    statemat[ii,ii] /= mean(sizes[uniques.==ii])
end
statemat3 = statemat ./ sum(statemat,dims=2)

open(datadir("transitionmatrix","$prefix-τ$τ-transitionmatrix-sizenormalized.txt"),"w") do io
    write(io,"$(statemat3)")
end

statemat = zeros(nstates,nstates)

for it = 1 : (length(uniques)-1)
    x = uniques[it]
    y = uniques[it+1]
    if !iszero(x) && !iszero(y); statemat[x,y] += 1 end
end

for it = 1 : (length(uniques)-2)
    x = uniques[it]
    y = uniques[it+1]
    z = uniques[it+2]
    if iszero(x) && iszero(z) && !iszero(y); statemat[y,y] += 1 end
end

statemat4 = statemat ./ sum(statemat,dims=2)

open(datadir("transitionmatrix","$prefix-τ$τ-transitionmatrix-uniques.txt"),"w") do io
    write(io,"$(statemat4)")
end