using DrWatson
@quickactivate "S2DTMatrices"

using DelimitedFiles
using ARMLive
using Statistics

τ = parse(Int64, ARGS[4])
prefix = "$(ARGS[1])-$(ARGS[2])-$(ARGS[3])-τ$τ"
fnc = "$prefix.nc"
ds  = NCDataset(datadir("states","$(ARGS[5])",fnc))
dt  = ds["time"][:] .- Hour(6)
states = ds["states"][:]; nstates = maximum(states)
uniques = ds["unique"][:]
sizes = ds["sizes"][:]
close(ds)

mkpath(datadir("timescales","$(ARGS[5])"))

count = zeros(nstates,nstates)
szmat = zeros(nstates,nstates)

for it = 1 : (length(uniques)-1)
    x = uniques[it]
    y = uniques[it+1]
    if !iszero(x) && !iszero(y)
        count[x,y] += 1
        szmat[x,y] += sizes[it]
    end
end

szmat = szmat ./ count

open(datadir("timescales","$(ARGS[5])","$prefix-raw.txt"),"w") do io
    write(io,"$(szmat)")
end