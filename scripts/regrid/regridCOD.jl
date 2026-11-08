using DrWatson
@quickactivate "S2DTMatrices"

using Printf
using GOESatellites
using NASAPrecipitation
using NCDatasets
using Statistics

yr = parse(Int,ARGS[1])
gds16 = GOESDataset("ABI-L2-COD",sector="C",satellite=16,path=datadir())
gds19 = GOESDataset("ABI-L2-COD",sector="C",satellite=19,path=datadir())
geo  = GeoRegion("SGP_LARGE",path=srcdir())
gvar = "COD"

minlon = floor(geo.W,digits=1) .- 0.05
maxlon =  ceil(geo.E,digits=1) .+ 0.05
minlat = floor(geo.S,digits=1) .- 0.05
maxlat =  ceil(geo.N,digits=1) .+ 0.05

dtvec = Date(yr) : Day(1) : Date(yr,12,31); ndays = length(dtvec)

ggrd = RegionGrid(gds19,geo); nlon,nlat = size(ggrd.lon)
rgln = zeros(Int,nlon,nlat)
rglt = zeros(Int,nlon,nlat)
rlon = minlon : 0.1 : maxlon; nrlon = length(rlon)
rlat = minlat : 0.1 : maxlat; nrlat = length(rlat)

for ilat = 1 : nlat, ilon = 1 : nlon
    rgln[ilon,ilat] = argmin(abs.(ggrd.lon[ilon,ilat] .- rlon))
    rglt[ilon,ilat] = argmin(abs.(ggrd.lat[ilon,ilat] .- rlat))
end

cellidx = [findall((rgln .== ilon) .& (rglt .== ilat)) for ilon = 1 : nrlon, ilat = 1 : nrlat]

rCOD = zeros(nrlon,nrlat,288) * NaN

for (ii,idt) in enumerate(dtvec)

    ds = idt < Date(2024,10,15) ? read(gds16,geo,"COD",idt,throw=false) : read(gds19,geo,"COD",idt,throw=false)
    rCOD .== NaN

    if !isnothing(ds)
        tCOD = nomissing(ds[gvar][:,:,:],0)
        dt   = ds["time"].var[:]; ndt = length(dt)
        close(ds)

        for it = 1 : ndt
            if !iszero(dt[it])
                iCOD = @view tCOD[:,:,it]
                for ilat = 1 : nrlat, ilon = 1 : nrlon
                    idx = cellidx[ilon,ilat]
                    if !isempty(idx)
                        iiCOD = view(iCOD, idx); iiCOD = iiCOD[.!isnan.(iiCOD)]
                        rCOD[ilon,ilat,it] = isempty(iiCOD) ? NaN : mean(iiCOD)
                    end
                end
            end
        end

        mkpath(datadir("regridded","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))"))
        fnc = datadir("regridded","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))","COD-$(geo.ID)-$(Dates.format(idt,dateformat"yyyymmdd")).nc")

        if isfile(fnc); rm(fnc,force=true) end
        ds = NCDataset(fnc,"c")

        defDim(ds,"longitude",nrlon)
        defDim(ds,"latitude",nrlat)
        defDim(ds,"time",288)

        nct  = defVar(ds,"time",Float64,("time",),attrib=Dict(
            "units"     => "minutes since 2000-01-01 00:00:00.0",
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

        nct.var[:] = dt
        nclon[:] = rlon
        nclat[:] = rlat
        ncCOD[:,:,:] = rCOD

        close(ds)

    end

end