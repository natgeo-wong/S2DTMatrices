using DrWatson
@quickactivate "S2DTMatrices"

using ARMLive
using Statistics

ads1 = ARMDataset(
    stream="sgpprecipbetseriesC1.c1",path=datadir(),
    start=Date(2011),stop=Date(2025,5,31)
)
ads2 = ARMDataset(
    stream="sgpprecipbetseriesC1.c0",path=datadir(),
    start=Date(2025,6),stop=Date(2025,12,31)
)

dtvec = ads1.start : Day(1) : ads2.stop; ndt = length(dtvec)
prcp  = zeros(144,ndt) * NaN

for (idt, dt) in enumerate(dtvec)
    @info "$(now()) - Reading data for $dt"; flush(stderr)
    ds1 = read(ads1,dt,throw=false)
    if !isnothing(ds1)
        dt = ds1["time"].var[:]
        pp = nomissing(ds1["precipitation_rate"][:],NaN)
        close(ds1)
        for it = 1 : 144
            iit = (dt .>= (it-1)*600) .& (dt .< it*600)
            iprcp = @views pp[iit]
            notnan = .!isnan.(iprcp)
            if sum(notnan) > 0
                prcp[it,idt] = mean(iprcp[notnan])
            end
        end
    else
        ds2 = read(ads2,dt,throw=false)
        if !isnothing(ds2)
            dt = ds2["time"].var[:]
            pp = nomissing(ds2["precipitation_rate"][:],NaN)
            close(ds2)
            for it = 1 : 144
                iit = (dt .>= (it-1)*600) .& (dt .< it*600)
                iprcp = @views pp[iit]
                notnan = .!isnan.(iprcp)
                if sum(notnan) > 0
                    prcp[it,idt] = mean(iprcp[notnan])
                end
            end
        end
    end
end

prcp[prcp.>100] .= NaN

fnc = datadir("RAIN-$(ads1.stream)-$(Dates.format(ads1.start,dateformat"yyyymmdd"))-$(Dates.format(ads2.stop,dateformat"yyyymmdd")).nc")
if isfile(fnc); rm(fnc,force=true) end
nds = NCDataset(fnc,"c")

defDim(nds,"time",ndt*144)

nct = defVar(nds,"time",Float64,("time",),attrib=Dict(
    "units"     => "minutes since $(ads1.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncprcp = defVar(nds,"precipitation_rate",Float64,("time",),attrib=Dict(
    "description"   => "Precipitation Rate",
    "standard_name" => "precipitation_rate",
    "units"         => "mm hr**-1"
))

nct.var[:] = collect(0.5 : (144*ndt)) * 10
ncprcp[:] = prcp
close(nds)