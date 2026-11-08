using DrWatson
@quickactivate "S2DTMatrices"

using DelimitedFiles
using ARMLive
using StatsBase

ads = ARMDataset(
    stream="sgpmfrsrcldod1minC1.c1",path=datadir(),
    start=Date(2011,1),stop=Date(2025,12,31)
)

dtvec = ads.start : Day(1) : ads.stop; ndays = length(dtvec)

COD = zeros(Float32,144,ndays) * NaN

for (ii,idt) in enumerate(dtvec)

    ds = read(ads,idt,throw=false)
    if !isnothing(ds)
        nt = ds.dim["time"]
        tCOD = nomissing(ds["optical_depth_instantaneous"][:],0)
        dt = ds["time"].var[:]
        close(ds)

        for it = 1 : 144
            iit = (dt .>= (it-1)*600) .& (dt .< it*600)
            iCOD = @views tCOD[iit]
            notnan = .!isnan.(iCOD)
            if sum(notnan) > 0
                COD[it,ii] = mean(iCOD[notnan])
            end
        end
    end

end

COD = COD[:]
dt = collect(0.5 : (144*ndays)) * 10; nt = length(dt)

fnc = datadir("COD-$(ads.stream)-$(Dates.format(ads.start,dateformat"yyyymmdd"))-$(Dates.format(ads.stop,dateformat"yyyymmdd")).nc")
if isfile(fnc); rm(fnc,force=true) end
ds = NCDataset(fnc,"c")

defDim(ds,"time",nt)

nct  = defVar(ds,"time",Float64,("time",),attrib=Dict(
    "units"     => "minutes since $(ads.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncCOD = defVar(ds,"COD",Float64,("time",),attrib=Dict(
    "description"   => "Cloud Optical Depth",
    "standard_name" => "cloud_optical_depth",
    "units"         => "dimensionless"
))

nct.var[:] = dt
ncCOD[:] = COD

close(ds)