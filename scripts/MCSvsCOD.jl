using DrWatson
@quickactivate "S2DTMatrices"

using Printf
using GOESatellites
using NASAPrecipitation
using NCDatasets
using StatsBase

dtvec = Date(2017,8,1) : Day(1) : Date(2025,12,31); ndays = length(dtvec)
geo = GeoRegion("SGP_LARGE",path=srcdir())

CODbin = -5 : 0.01 : 3
COD1 = zeros(Int,length(CODbin)-1)
COD2 = zeros(Int,length(CODbin)-1)
COD3 = zeros(Int,length(CODbin)-1)
COD4 = zeros(Int,length(CODbin)-1)

for (ii,idt) in enumerate(dtvec)

    fncCOD = datadir("regridded-hourly","$(geo.ID)","$(year(idt))","$(@sprintf("%02d",month(idt)))","COD-$(geo.ID)-$(Dates.format(idt,dateformat"yyyymmdd")).nc")
    fncMCS = datadir("MCStracks","fengzhe21","$(geo.ID)","$(year(idt))","$(Dates.format(idt,dateformat"yyyymmdd")).nc")

    if isfile(fncCOD) && isfile(fncMCS)
        
        dsMCS = NCDataset(fncMCS)
        cldtyp = dsMCS["cloudtype"][:,:,:]
        close(dsMCS)
        
        dsCOD = NCDataset(fncCOD)
        COD = log10.(dsCOD["COD"][2:(end-1),2:85,:])
        close(dsCOD)

        h = fit(Histogram, (@view COD[cldtyp .== 1]), CODbin); COD1[:] += h.weights
        h = fit(Histogram, (@view COD[cldtyp .== 2]), CODbin); COD2[:] += h.weights
        h = fit(Histogram, (@view COD[cldtyp .== 3]), CODbin); COD3[:] += h.weights
        h = fit(Histogram, (@view COD[cldtyp .== 4]), CODbin); COD4[:] += h.weights

    end

end

open(datadir("MCSvsCOD1.txt"),"w") do io
    write(io,"$(COD1)")
end
open(datadir("MCSvsCOD2.txt"),"w") do io
    write(io,"$(COD2)")
end
open(datadir("MCSvsCOD3.txt"),"w") do io
    write(io,"$(COD3)")
end
open(datadir("MCSvsCOD4.txt"),"w") do io
    write(io,"$(COD4)")
end