using Images
using ImageSegmentation
using ImageFiltering
using CSV
using DataFrames, DataFramesMeta
using Dates
using Statistics, StatsBase, Random
using Plots, StatsPlots, Plots.PlotMeasures

function to_array(echodf)
    echo_arr = unstack(echodf, :Interval, :Layer_depth_min, :Sv_mean)
    echo_arr = Matrix(echo_arr[:, 2:end])
    echo_arr[ismissing.(echo_arr) .| (echo_arr .< -90)] .= -90
    return disallowmissing(echo_arr)
end

function find_blobs(echo, threshold)
    bw = echo .> threshold
    markers = label_components(bw)
    w = watershed(echo, markers)
    return labels_map(w)
end

function segment_stats(label, labelmap, echo)
    indices = findall(x -> x == label, labelmap)
    bottom, top = extrema(i.I[2] for i in indices)
    left, right = extrema(i.I[1] for i in indices)
    height = top - bottom + 1
    @assert height > 0 println("bottom: $bottom, top: $top")
    width = right - left + 1
    values = [echo[i] for i in indices]
    npixels = length(values)
    μ = mean(values)
    σ = std(values)
    return (bottom=bottom, top=top, left=left, right=right,
        height=height, width=width, npixels=npixels,
        mean=μ, std=σ, total=μ*npixels)
end

function process_day(df)
    echo_arr = to_array(df)
    datetimes = sort(unique(df.datetime))
    depths = sort(unique(df.Layer_depth_min))
    markers = find_blobs(echo_arr, -60)
    schools = DataFrame(segment_stats(i, markers, exp10.(echo_arr/10))
        for i in unique(markers))
    schools[!, :school_id] .= unique(markers)
    schools = @transform(schools,
        # school_id .= unique(markers),
        :start_time = datetimes[:left],
        # top and bottom are reversed to go from image to depth coordinates
        :bottom = depths[:top],
        :top = depths[:bottom])
    schools = @transform(schools,
        :middle = (:bottom .+ :top) ./ 2,
        :date = fill(first(df.date), nrow(schools)))

    schools = @subset(schools, 2 .< :npixels .< 1000, :width .<= 12)

    # now calculate center of mass, ignoring regions classified as schools
    sv = exp10.(echo_arr/10)
    ii = in(schools.school_id).(markers)
    sv[ii] .= 0
    cm = vec(sum(sv .* depths', dims=2) ./ sum(sv, dims=2))
    layertop = depths[vec(mapslices(ping -> findfirst(cumsum(ping) .> 1e-6), sv, dims=2))]
    layerbase = reverse(depths)[vec(mapslices(ping -> findfirst(cumsum(reverse(ping)) .> 1e-6), sv, dims=2))]
    center_of_mass = DataFrame(date=first(df.date), datetime=datetimes,
        cm=cm, layertop=layertop, layerbase=layerbase)
    return schools, center_of_mass
end

function get_random_color(seed)
    Random.seed!(seed)
    rand(RGB{N0f8})
end

#######################################################################################
# This section illustrates the school-detection procedure on a three-day subset of data.
#######################################################################################
# Reading in a month of acoustic data and formatting it
echo = CSV.read(joinpath(@__DIR__, "../data/monthly/deimos-2019-11.csv"),
    DataFrame, normalizenames=true)

echo = @chain echo begin
    @transform(:datetime = Date.(string.(:Date_M), dateformat"yyyymmdd") .+ :Time_M .- Hour(8))
    @transform(:hour = hour.(:datetime),
              :date = Date.(:datetime))
end

# Selecting three days with some schools
echo_sub = @subset(echo,
    :datetime .> DateTime(2019, 11, 5, 0),
    :datetime .< DateTime(2019, 11, 8, 0),
    :Layer_depth_min .<300)
# convert the data frame to a 2d array
echo_arr = to_array(echo_sub)
datetimes = sort(unique(echo_sub.datetime))
depths = sort(unique(echo_sub.Layer_depth_min))

# locate bright regions on the echogram and label them
markers = find_blobs(echo_arr, -60)
heatmap(map(i->get_random_color(i), markers'))

# calculate statistics for each individual region, storing them in a DataFrame
schools = DataFrame(segment_stats(i, markers, exp10.(echo_arr/10))
  for i in unique(markers))
schools.school_id = unique(markers)

schools = @transform(schools,
  :start_time = datetimes[:left],
  # top and bottom are reversed to go from image to depth coordinates
  :bottom = depths[:top],
  :top = depths[:bottom])
schools = @transform(schools,
  :middle = (:bottom .+ :top) ./ 2,
  :date = fill(first(echo_sub.date), nrow(schools)))

# Drop schools that are too small or too big
schools = @subset(schools, 2 .< :npixels .< 1000, :width .< 12)

# Blank out schools from the echogram
sv = exp10.(echo_arr/10)
ii = in(schools.school_id).(markers)
sv[ii] .= 0
# Echogram with the schools blanked out
heatmap(10log10.(sv'), yflip=true)

# Calculate time series of several layer metrics, ignoring regions classified as schools
cm = vec(sum(sv .* depths', dims=2) ./ sum(sv, dims=2))
center_of_mass = DataFrame(date=first(echo_sub.date), datetime=datetimes, cm=cm)
layertop = depths[vec(mapslices(ping -> findfirst(cumsum(ping) .> 1e-6), sv, dims=2))]
layerbase = reverse(depths)[vec(mapslices(ping -> findfirst(cumsum(reverse(ping)) .> 1e-6), sv, dims=2))]

heatmap(1:432, depths, 10log10.(sv)', yflip=true)
plot!(cm, color=:white, label = "")
plot!(layertop, color=:green, label="")
plot!(layerbase, color=:blue, label="")

## Process all DEIMOS data
schools_list = []
cm_list = []
files = readdir(joinpath(@__DIR__, "../data/monthly"), join=true)

for f in files
    println("Reading $(basename(f))...")
    echo = CSV.read(f, DataFrame, normalizenames=true)
    echo = @chain echo begin
        @subset(15 .<= :Layer_depth_min .<= 500)
        @transform(:datetime = Date.(string.(:Date_M), dateformat"yyyymmdd") .+ :Time_M .- Hour(8))
        @transform(:hour = hour.(:datetime),
                  :date = Date.(:datetime))
    end
    println("Detecting schools...")
    schools = []
    center_of_mass = []
    for echo_day in groupby(echo, :date)
        schools_day, center_of_mass_day = process_day(echo_day)
        push!(schools, schools_day)
        push!(center_of_mass, center_of_mass_day)
    end
    schools = vcat(schools...)
    schools[!, :school_id] .= string.(schools.date) .* "-S" .* string.(schools.school_id)
    center_of_mass = vcat(center_of_mass...)

    push!(schools_list, schools)
    push!(cm_list, center_of_mass)
end
schools = vcat(schools_list...)
center_of_mass = @chain vcat(cm_list...) begin
    @transform(:hour = hour.(:datetime))
end

histogram(log10.(schools.total))
histogram(schools.bottom, bins=200)
histogram(@subset(schools, :width .<= 2).bottom, bins=100)
schools = @subset(schools, :bottom .< 200)

schools_summary = @chain schools begin
    groupby(:date)
    @combine(:nschools = length(:school_id),
        :sa = sum(:total),
        :mean_depth = mean(:middle),
        :max_depth = maximum(:bottom))
end

alldates = DataFrame(date=first(schools.date):Day(1):last(schools.date))
schools_summary = leftjoin(alldates, schools_summary, on=:date)
schools_summary = @chain schools_summary begin
    @transform(:nschools = replace(:nschools, missing => 0),
               :sa = replace(:sa, missing => 0))
    @transform(:Sa = 10log10.(:sa))
    @orderby(:date)
end

@df schools_summary plot(:date, :sa)
@df schools_summary plot(:date, :nschools)
@df schools_summary plot(:date, [:mean_depth, :max_depth])
@df schools_summary scatter(:nschools, :Sa)

@df center_of_mass plot(:datetime, :cm)

center_of_mass = @chain center_of_mass begin
    @transform(:cm_smooth = mapwindow(mean, :cm, 6*3+1),
               :layertop_smooth = mapwindow(mean, :layertop, 6*3+1),
               :layerbase_smooth = mapwindow(mean, :layerbase, 6*3+1))
end

@df center_of_mass[10_000:11_000,:] plot(:datetime,
    [:cm, :cm_smooth, :layertop, :layertop_smooth, :layerbase, :layerbase_smooth])

CSV.write(joinpath(@__DIR__, "..\\data\\schools.csv"), schools)
CSV.write(joinpath(@__DIR__, "..\\data\\schools_summary.csv"), schools_summary)
CSV.write(joinpath(@__DIR__, "..\\data\\center_of_mass.csv"), center_of_mass)

schools = CSV.read(joinpath(@__DIR__, "..\\data\\schools.csv"), DataFrame)
schools_summary = CSV.read(joinpath(@__DIR__, "..\\data\\schools_summary.csv"), DataFrame)
center_of_mass = CSV.read(joinpath(@__DIR__, "..\\data\\center_of_mass.csv"), DataFrame)


## Calculating migration statistics

noon(hour) = hour in [12, 13, 14]
midnight(hour) = hour in [0, 1, 2]

function hourlabel(hour)
    if noon(hour)
        return "noon"
    elseif midnight(hour)
        return "midnight"
    else
        return "other"
    end
end

migration = @chain center_of_mass begin
    @transform(:hourlabel = hourlabel.(:hour))
    @subset(:hourlabel .!= "other")
    groupby([:date, :hourlabel])
    @combine(:cm=median(:cm_smooth),
        :layertop=median(:layertop_smooth))
end

migration = innerjoin(unstack(migration, :date, :hourlabel, :cm, renamecols=x->Symbol(:cm_, x)),
     unstack(migration, :date, :hourlabel, :layertop, renamecols=x->Symbol(:layertop_, x)),
     on=:date)
migration = @chain migration begin
    @transform(:Δcm = :cm_noon .- :cm_midnight)
    @transform(:Δlayertop = :layertop_noon .- :layertop_midnight)
    leftjoin(schools_summary, on=:date)
    @transform(:nschools = replace(:nschools, missing => 0),
               :sa = replace(:sa, missing => 0),
               :Sa = replace(:Sa, -Inf => missing))
end

@df migration plot(:date, [:Δcm, :Δlayertop], labels=["ΔCM" "ΔLayertop"])

@df @subset(center_of_mass, Date(2019, 4, 1) .< :date .< Date(2019, 5, 31)) plot(:datetime, :cm_smooth)

plot(plot(migration.date, migration.layertop_noon, ylabel="Migration"),
    plot(migration.date, migration.sa, ylabel="schools"), layout=(2,1))

@df migration scatter(:nschools, [:cm_noon, :layertop_noon], labels=["CM" "Layer top"])
@df migration scatter(:Sa, [:cm_noon, :layertop_noon], labels=["CM" "Layer top"])

@df migration scatter(:nschools, [:Δcm, :Δlayertop], labels=["CM" "Layer top"])
@df migration scatter(:Sa, [:Δcm, :Δlayertop], labels=["CM" "Layer top"])


## Regression modeling
using GLM

# CMigration, center of mass
m1 = lm(@formula(Δcm ~ nschools), migration)
r2(m1)
m2 = lm(@formula(Δcm ~ Sa), migration)
r2(m2)
m3 = lm(@formula(Δcm ~ Sa + nschools), migration)
r2(m3)

# Migration, layer top
m1 = lm(@formula(Δlayertop ~ nschools), migration)
r2(m1)
m2 = lm(@formula(Δlayertop ~ Sa), migration)
r2(m2)
m3 = lm(@formula(Δlayertop ~ Sa + nschools), migration)
r2(m3)


p1 = @df migration scatter(:Sa, :Δlayertop, legend=false, markerstrokewidth=0, color=:grey,
    xlabel="School backscatter (dB re m² minute)", ylabel="Migration magnitude (m)",
    title="A", dpi=600)
plot!(p1, Sa -> coef(m2)[1] + coef(m2)[2] * Sa, -58, -12, color=:black)

# Noon CM depth
m1 = lm(@formula(cm_noon ~ nschools), migration)
r2(m1)
m2 = lm(@formula(cm_noon ~ Sa), migration)
r2(m2)
m3 = lm(@formula(cm_noon ~ Sa + nschools), migration)
r2(m3)

# Noon layer top depth
m1 = lm(@formula(layertop_noon ~ nschools), migration)
r2(m1)
m2 = lm(@formula(layertop_noon ~ Sa), migration)
r2(m2)
m3 = lm(@formula(layertop_noon ~ Sa + nschools), migration)
r2(m3)

p2 = @df migration scatter(:Sa, :layertop_noon, legend=false, yflip=true,
    markerstrokewidth=0, color=:grey, title="B",
    xlabel="School backscatter (dB re m² minute)", ylabel="Depth at noon (m)", dpi=600)
plot!(p2, Sa -> coef(m2)[1] + coef(m2)[2] * Sa, -58, -12, color=:black)

p = plot(p1, p2, size=(1000, 400), dpi=600, title_loc=:left, margin=5mm)
savefig(p, joinpath(@__DIR__, "../graphics/schools_vs_layers.png"))

## Check against oceanography too

oceanography = CSV.read(joinpath(@__DIR__, "../data/oceanography.csv"), DataFrame)
oceanography = @chain oceanography begin
    @subset(.! ismissing.(:date))
    groupby(:date)
    @combine(:temperature = mean(skipmissing(:temperature)),
             :upwelling = mean(skipmissing(:upwelling)),
             :sealevel = mean(skipmissing(:sealevel_lowpass)))
    @transform(upwelling = replace(:upwelling, NaN => missing))
    rightjoin(migration, on=:date)
endlibrary(lubridate)
library(viridis)
library(scales)

@df oceanography scatter(:temperature, :layertop_noon, legend=false,
    xlabel="SST (°C)", ylabel="Noon depth (m)", dpi=300)
savefig(joinpath(@__DIR__, "../graphics/temperature_vs_layertop.png"))
@df oceanography scatter(:sealevel, :layertop_noon, legend=false,
    xlabel="Sea level (m)", ylabel="Noon depth (m)", dpi=300)
savefig(joinpath(@__DIR__, "../graphics/sealevel_vs_layertop.png"))
@df oceanography scatter(:upwelling, :layertop_noon, legend=false,
    xlabel="Upwelling (kg m⁻¹ s⁻¹)", ylabel="Noon depth (m)", dpi=300)
savefig(joinpath(@__DIR__, "../graphics/upwelling_vs_layertop.png"))

mtemp = lm(@formula(layertop_noon ~ temperature), oceanography)
r2(mtemp)
msealevel = lm(@formula(layertop_noon ~ sealevel), oceanography)
r2(msealevel)
mupwelling = lm(@formula(layertop_noon ~ upwelling), dropmissing(oceanography))
r2(mupwelling)
