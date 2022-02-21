#using ArgCheck
using Dates, SparseArrays, StatsBase
using CSV
using Pandas, DataFrames#master
using Dates
using HTTP
#using PooledArrays
using Tables
using TypedTables
using AxisKeys
using Parameters
#using Optim
using LinearAlgebra, Distributions
using ArraysOfArrays
DataFrame = DataFrames.DataFrame

in_days(x::Day) = x.value
in_days(x::TimePeriod) = Day(x).value
export in_days


# TODO: Make sure all time-dependent data starts at least on September 1st of 2020!

function Base.:*(::Float64, ::String) return NaN end
function Base.:+(::String, ::Float64) return NaN end

function download_file(url::AbstractString, filename::AbstractString)
    if !isfile(filename)
        @info "Downloading \"$filename\" from \"$url\""
        open(filename, "w") do io
            content = HTTP.get(url).body
            content_cmp = HTTP.get(url).body
            if (content == content_cmp)
                write(io, content)
                nothing
            else
                throw(ErrorException("Invalid/broken download"))
            end
        end
    else
        # @info "File \"$filename\" already present, keeping it."
    end
    return filename
end
export download_file


getcol(df::DataFrame, colname::AbstractString) = getproperty(df, Symbol(colname))

getcol(f::Base.Callable, df::DataFrame, colname::AbstractString) = passmissing(f).(getcol(df, colname))



# ToDo: Use CodecZlib to gzip data before writing:
rki_data_path() = download_file("https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data", joinpath("data", "RKI_COVID19.csv"))

"""
    rki_data()::TypedTables.Table

[RKI Covid-19 data](https://hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6)
data.
"""
function rki_data()
    time_expr1 = r"^([0-9][0-9][0-9][0-9])/([0-9][0-9])/([0-9][0-9]) ([0-9][0-9]):([0-9][0-9]):([0-9][0-9])$"
    time_expr2 = r"^([0-9][0-9])\.([0-9][0-9])\.([0-9][0-9][0-9][0-9]), ([0-9][0-9]):([0-9][0-9]) Uhr$"
    
    function conv_time(s::AbstractString)
        if s[5] == '/'
            Y, m, d, H, M, S = parse.(Int, match(time_expr1, s).captures)
            # DateTime(Y, m, d, H, M, S)
            Date(Y, m, d)
        else
            d, m, Y, H, M = parse.(Int, match(time_expr2, s).captures)
            # DateTime(Y, m, d, H, M)
            Date(Y, m, d)
        end
    end
    
    df = CSV.File(rki_data_path()) |> DataFrame
    # show(df, allcols = true)
    df.Meldedatum = conv_time.(df.Meldedatum)
    df.Datenstand = conv_time.(df.Datenstand)
    df.Refdatum = conv_time.(df.Refdatum)
    df.IstErkrankungsbeginn = Bool.(df.IstErkrankungsbeginn)
    ;

    # show(propertynames(df))

    TypedTables.Table(
        fid = df.FID::Vector{Int64},
        report_date = df.Meldedatum::Vector{Dates.Date},
        ref_date = df.Refdatum::Vector{Dates.Date},
        is_onset_date = df.IstErkrankungsbeginn::BitVector,
        state_id = df.IdBundesland::Vector{Int64},
        state = String.(df.Bundesland),
        district_id = df.IdLandkreis::Vector{Int64},
        district = String.(df.Landkreis),
        agegroup = String.(df.Altersgruppe),
        gender = String.(df.Geschlecht),
        n_cases = df.AnzahlFall::Vector{Int64},
        n_deaths = df.AnzahlTodesfall::Vector{Int64},
        n_healed = df.AnzahlGenesen::Vector{Int64},
        new_case = df.NeuerFall::Vector{Int64},
        new_death = df.NeuerTodesfall::Vector{Int64},
        new_healed = df.NeuGenesen::Vector{Int64},
        data_date = df.Datenstand::Vector{Dates.Date},
    )
end
export rki_data

age_groups() = (0:8)*10
export age_groups

regions_ger() = ["BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV", "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH"]
export regions_ger

function rki_data_proxy(rki;by_ref_date::Bool=false,isonset::Bool=false,exclusive::Bool=false)
    if isonset
        rki = rki[findall(rki.is_onset_date)]
    elseif exclusive
        rki = rki[findall(Vector{Bool}((rki.is_onset_date .- 1 ) * (-1)))]
    end

    dates = rki.report_date # rki.to Meldedatum; rki.ref_date Referenzdatum; rki.is_from::bool IST Erkrankungsbeginn
    if by_ref_date
        dates = rki.ref_date
    end
    global_from, global_until = sort(unique(dates))[[1,end]]
    global_dates = global_from:Day(1):global_until
#     dates_id = Vector{Int}( (dates - global_from ) / Day(1) .+ 1)
    global_genders = sort(unique(rki.gender))
    gender_id = (rki.gender .== "M") .+1 #TODO ignores existence of type unknown
    global_ages = [ 0, 5, 15, 35, 60,  80, -1]
    global_agegroups = sort(unique(rki.agegroup))
    agegroup_id = zeros(Int,length(rki.agegroup))
    for i in 1:length(global_agegroups)
        agegroup_id += i .* (rki.agegroup .== global_agegroups[i])
    end

    bl = ["SH", "HH", "NI", "HB", "NW", "HE", "RP", "BW", "BY", "SL", "BE", "BB", "MV",  "SN", "ST", "TH" ]

    global_regions = unique(sort(rki.state_id))
    dims = Tuple([length(global_dates), length(global_regions), length(global_agegroups)])
    cases = KeyedArray(zeros(Int, dims...), date=Vector(global_dates), bl=bl, age=global_ages)
    deaths = copy(cases)
    healed = copy(cases)

    for i in 1:length(dates)
        a = agegroup_id[i]
        r = rki.state_id[i]
        cases(dates[i])[r,a] += rki.n_cases[i]
        deaths(dates[i])[r,a] += rki.n_deaths[i]
        healed(dates[i])[r,a] += rki.n_healed[i]
    end
    return ( cases = copy(cases(regions_ger())), deaths = copy(deaths(regions_ger())), healed = copy(healed(regions_ger())) )
end

function rki_data_proxy(by_ref_date::Bool=false,isonset::Bool=false,exclusive::Bool=false)
    rki = rki_data()
    rki_data_proxy(rki,by_ref_date,isonset,exclusive)
end
export rki_data_proxy



"""
    download_all_data()

Download all required data to the local directory "data".
"""
function download_all_data()
    # ToDo: Add option to update already downloaded data.
    
    @sync begin

        # National data:
        @async rki_data_path()

    end

    return nothing
end
export download_all_data


