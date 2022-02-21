using Pkg
Pkg.activate(abspath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using ShowMe
using BAT4Epidemiology
using BAT
using Distributions
using Plots
using AxisKeys

gr(size=(1.1*850, 1.1*600), thickness_scaling = 1.3)
plot = Plots.plot

include("../src/CovidData.jl")
plot = Plots.plot
#===============================================================#

#from  = Date(2020,10,1); until = Date(2021, 7,1);
#t_rng = 0.0:210;
#dates = from:Day(1):until

#rki = rki_data();
#@unpack cases, deaths = rki_data_proxy(rki;by_ref_date=true,isonset=true,exclusive=true);
#cases_onset = copy(cases(dates));
#@unpack cases, deaths = rki_data_proxy(rki;by_ref_date=false,isonset=true,exclusive=true);
#cases_delayed_onset = copy(cases(dates));

