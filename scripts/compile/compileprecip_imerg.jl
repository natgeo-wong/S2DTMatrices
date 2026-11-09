using DrWatson
@quickactivate "S2DTMatrices"

using DelimitedFiles, Statistics
using NASAPrecipitation, RegionGrids
using StatsBase

sname = "SGP"; sdata = readdlm(srcdir("armstations.txt"),',')
sID = sdata[:,1]
slon = sdata[findfirst(sID .== sname),2]
slat = sdata[findfirst(sID .== sname),3]

npd = IMERGFinalHH(start=Date(2011),stop=Date(2025,12,31),path=datadir())
geo = GeoRegion("$(sname)_2x2",path=srcdir())
lsd = getLandSea(npd,geo)

ggrd = RegionGrid(geo,lsd.lon,lsd.lat)

ind1 = nearest(Point2(slon,slat),ggrd)
ind4 = nearest(Point2(slon,slat),ggrd,n=4); d4 = distances(Point2(slon,slat),ggrd,n=4)
wgts = weights((1 ./d4).^2)
dall = distances(Point2(slon,slat),ggrd)
iall = dall .< (100e3/sqrt(pi))

dtvec = npd.start : Day(1) : npd.stop; ndt = length(dtvec)
prcp  = zeros(48,ndt) * NaN
prcp1 = zeros(48,ndt) * NaN
prcp4 = zeros(48,ndt) * NaN
p100k  = zeros(48,ndt) * NaN
prcpt = zeros(Float32,length(lsd.lon),length(lsd.lat),48)

for (idt,dt) in enumerate(dtvec)

    ds = read(npd,geo,dt,throw=false)
    if !isnothing(ds)
        NCDatasets.load!(ds["precipitation"].var,prcpt,:,:,:)
        close(ds)
        for it = 1 : 48
            iprcp = @views prcpt[:,:,it]
            prcp1[it,idt] = iprcp[ind1]
            prcp4[it,idt] = mean(iprcp[ind4],wgts)
            p100k[it,idt] = mean(iprcp[iall])
        end
    end

end

fnc = datadir("RAIN-$(npd.ID)-$(geo.ID)-$(Dates.format(npd.start,dateformat"yyyymmdd"))-$(Dates.format(npd.stop,dateformat"yyyymmdd")).nc")

if isfile(fnc); rm(fnc,force=true) end
nds = NCDataset(fnc,"c")

defDim(nds,"time",ndt*24)

nct = defVar(nds,"time",Float64,("time",),attrib=Dict(
    "units"     => "minutes since $(npd.start) 00:00:00.0",
    "long_name" => "time",
    "calendar"  => "gregorian",
))

ncprcp1 = defVar(nds,"prcp_nn",Float64,("time",),attrib=Dict(
    "description"   => "Precipitation (Nearest Neighbour)",
    "units"         => "mm"
))

ncprcp2 = defVar(nds,"prcp_4",Float64,("time",),attrib=Dict(
    "description"   => "Precipitation (Weighted-Mean (1/r^2) 4 Nearest Neighbours)",
    "units"         => "mm"
))

ncprcp3 = defVar(nds,"prcp_100km2",Float64,("time",),attrib=Dict(
    "description"   => "Precipitation (Mean 100km2 Circle)",
    "units"         => "mm"
))

nct[:] = collect(0.5:(ndt*24)) * 30
ncprcp1[:] = prcp1[:] * 3600
ncprcp2[:] = prcp4[:] * 3600
ncprcp3[:] = p100k[:] * 3600
close(nds)