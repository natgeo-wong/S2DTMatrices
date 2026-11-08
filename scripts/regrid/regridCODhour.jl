using DrWatson
@quickactivate "S2DTMatrices"

using Printf
using GOESatellites
using NASAPrecipitation
using NCDatasets
using Statistics

yr = parse(Int,ARGS[1])
geo  = GeoRegion("SGP_LARGE",path=srcdir())
gvar = "COD"

minlon = floor(geo.W,digits=1) .- 0.05
maxlon =  ceil(geo.E,digits=1) .+ 0.05
minlat = floor(geo.S,digits=1) .- 0.05
maxlat =  ceil(geo.N,digits=1) .+ 0.05

dtvec = Date(yr) : Day(1) : Date(yr,12,31); ndays = length(dtvec)

rlon = minlon : 0.1 : maxlon; nlon = length(rlon)
rlat = minlat : 0.1 : maxlat; nlat = length(rlat)

rCOD = zeros(nlon,nlat,24) * NaN

for (ii,idt) in enumerate(DateTime.(dtvec))
    rCOD .== NaN

    fnc = datadir("regridded","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))","COD-$(geo.ID)-$(Dates.format(idt,dateformat"yyyymmdd")).nc")

    if isfile(fnc)
        ds = NCDataset(fnc)
        tCOD = nomissing(ds[gvar][:,:,:],0)
        dt   = ds["time"][:]; ndt = length(dt)
        close(ds)

        for it = 1 : 24
            ii = (dt .>= (idt.+Hour(1)*(it-1))) .& (dt .< (idt.+Hour(1)*it))
            if any(ii)
                for ilat = 1 : nlat, ilon = 1 : nlon
                    iCOD = @view tCOD[ilon,ilat,ii]
                    notnan = .!isnan.(iCOD)
                    rCOD[ilon,ilat,it] = !iszero(sum(notnan)) ? mean(@view iCOD[notnan]) : NaN
                end
            end
        end

        mkpath(datadir("regridded-hourly","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))"))
        fnc = datadir("regridded-hourly","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))","COD-$(geo.ID)-$(Dates.format(idt,dateformat"yyyymmdd")).nc")

        if isfile(fnc); rm(fnc,force=true) end
        ds = NCDataset(fnc,"c")

        defDim(ds,"longitude",nrlon)
        defDim(ds,"latitude",nrlat)
        defDim(ds,"time",24)

        nct  = defVar(ds,"time",Float64,("time",),attrib=Dict(
            "units"     => "hours since $idt 00:00:00.0",
            "long_name" => "time",
            "calendar"  => "gregorian",
        ))

        nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
            "units"         => "degrees_east",
            "standard_name" => "longitude",
        ))

        nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
            "units"         => "degrees_north",
            "standard_name" => "latitude",
        ))

        ncCOD  = defVar(ds,"COD",Float64,("longitude","latitude","time",),attrib=Dict(
            "description"   => "Cloud Optical Depth",
            "standard_name" => "cloud_optical_depth",
            "units"         => "dimensionless"
        ))

        nct.var[:] = 0.5:24
        nclon[:] = rlon
        nclat[:] = rlat
        ncCOD[:,:,:] = rCOD

        close(ds)

    end

end