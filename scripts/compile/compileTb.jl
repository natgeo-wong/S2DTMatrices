using DrWatson
@quickactivate "S2DTMatrices"

using Dates, Logging, Printf, Statistics
using DataInterpolations
using NCDatasets
using GeoRegions, RegionGrids

geo  = GeoRegion("SGP_LARGE",path=srcdir())

minlon = floor(geo.W,digits=1) .+ 0.05
maxlon =  ceil(geo.E,digits=1) .- 0.05
minlat = floor(geo.S,digits=1) .+ 0.05
maxlat =  ceil(geo.N,digits=1) .- 0.05

dtvec = Date(2011) : Day(1) : Date(2025,12,31); ndays = length(dtvec)

rlon = minlon : 0.1 : maxlon; nlon = length(rlon)
rlat = minlat : 0.1 : maxlat; nlat = length(rlat)

ggrd = RegionGrid(geo,rlon,rlat)
ilon = argmin(abs.(rlon .- (geo.W + geo.E)/2))
ilat = argmin(abs.(rlat .- (geo.S + geo.N)/2))

ttb = zeros(24,ndays)*NaN

for (ii,idt) in enumerate(DateTime.(dtvec))
    
    yr = year(idt)
    dtstr = Dates.format(idt,dateformat"yyyymmdd")
    fnc = datadir("MCStracks","fengzhe21","$(geo.ID)","$(yr)","$(dtstr).nc")
    if isfile(fnc)
        
        @info "Reading Tb data for $idt"; flush(stderr)
        ds = NCDataset(fnc)
        ttb[:,ii] .= ds["tb"][ilon,ilat,:]
        close(ds)
        
    else

        @info "No Tb data found for $idt"; flush(stderr)

    end

end

tt = collect(30 : 60 : (ndays*60*24))
nt = collect( 5 : 10 : (ndays*60*24))
ttb = ttb[:]
ttb[isnan.(ttb)] .= 0

@info "Makima Interpolation"; flush(stderr)
aki = AkimaInterpolation(ttb,tt,modified=true,extrapolation=ExtrapolationType.Constant)
tb = aki.(nt)

fnc = datadir("TB-SGP-$(Dates.format(dtvec[1],dateformat"yyyymmdd"))-$(Dates.format(dtvec[2],dateformat"yyyymmdd")).nc")
if isfile(fnc); rm(fnc,force=true) end
ds = NCDataset(fnc,"c")

defDim(ds,"time",length(nt))

nct  = defVar(ds,"time",Float64,("time",),attrib=Dict(
    "units"     => "minutes since 2011-01-01 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncTB  = defVar(ds,"Tb",Float64,("time",),attrib=Dict(
    "description"   => "Brightness Temperature nearest to ARM Facility",
    "standard_name" => "brightness_temperature",
    "units"         => "K"
))

nct.var[:] = nt
ncTB[:] = tb

close(ds)