using DrWatson
@quickactivate "S2DTMatrices"

using Dates
using DelimitedFiles
using ARMLive

fnc = "COD-$(ARGS[1])-$(ARGS[2])-$(ARGS[3]).nc"
ds  = NCDataset(datadir(fnc))
dt  = ds["time"][:]; hr = hour.(dt .- Hour(6)); mo = month.(dt .- Hour(6))
COD = ds["COD"][:]; nCOD = length(COD)
close(ds)

fnc = "RAIN-imergv7finalhh-SGP_2x2-$(ARGS[2])-$(ARGS[3]).nc"
ds  = NCDataset(datadir(fnc))
tprcp = ds["prcp_100km2"][:]
close(ds)

precip = zeros(3,length(tprcp))
precip[1,:] .= tprcp
precip[2,:] .= tprcp
precip[3,:] .= tprcp
precip = reshape(precip,:)

τ = parse(Int64, ARGS[4])

states = zeros(nCOD)
states2 = zeros(nCOD) * NaN
sizes = zeros(nCOD)

states[iszero.(COD)] .= 1
states[(COD.>0).&(COD.<=τ)] .= 2
states[(COD.>τ).&(precip.>0)] .= 3
states[(COD.>τ).&iszero.(precip)] .= 4
states[(COD.>τ).&isnan.(precip)] .= 5
states[(hr .< 6) .| (hr .> 22)] .= 0
states[(mo .< 4) .| (mo .> 9)] .= 0

states = reshape(states,144,:); ndays = size(states,2)
for iday = 1 : ndays
    istates = @view states[:,iday]; i3 = istates .== 3
    first3 = findfirst(i3); last3 = findlast(i3)
    if !isnothing(first3) && !isnothing(last3)
        if last3 !== first3
            istates[first3:last3] .= 3
        end
    end
end
states = reshape(states,:)

global ii = 1
for it = 1 : nCOD
    states2[ii] = states[it]
    sizes[ii] += 1
    if it < nCOD
        if states[it+1] !== states[it]
            global ii += 1
        end
    end
end

states2 = states2[1:ii]
sizes = sizes[1:ii]

mkpath(datadir("states","f3l3"))
fnc = datadir("states","f3l3","COD-$(ARGS[1])-RAIN-imergv7finalhh-100km2-$(ARGS[2])-$(ARGS[3])-τ$τ-states.nc")
if isfile(fnc); rm(fnc,force=true) end
ds = NCDataset(fnc,"c")

defDim(ds,"time",nCOD)
defDim(ds,"states",ii)

nct  = defVar(ds,"time",Float64,("time",),attrib=Dict(
    "units"     => "minutes since 2000-01-01 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncstates = defVar(ds,"states",Int64,("time",),attrib=Dict(
    "description"   => "Cloud States",
    "standard_name" => "cloud_states",
    "units"         => "0=Missing,1=Clear,2=Shallow,3=Deep"
))

ncstates2 = defVar(ds,"unique",Int64,("states",),attrib=Dict(
    "description"   => "Unique Cloud States",
    "standard_name" => "unique_cloud_states",
    "units"         => "0=Missing,1=Clear,2=Shallow,3=Deep"
))

ncsizes = defVar(ds,"sizes",Int64,("states",),attrib=Dict(
    "description"   => "State Sizes",
    "standard_name" => "State sizes",
    "units"         => "minutes"
))

nct[:] = dt
ncstates[:] = states
ncstates2[:] = states2
ncsizes[:] = sizes * 10

close(ds)