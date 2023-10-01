# using Pkg
# Pkg.activate(".")
using Plots
using CSV, DataFrames, XLSX, StatsBase, MixedModels, Dates, Bootstrap
df = DataFrame(XLSX.readtable("Joint Microadenoma Project Data Collection PRLACTINOMAS.xlsx", 1))
# df = DataFrame(XLSX.readtable("microadenoma_all.xlsx", 1))


intervention = df[!, Symbol("Any intervention during follow-up? (0=None, 1=Surgery, 2=Medical therapy, 3=Radiation therapy, 4=other)")]
# ividx = intervention .== 0

prolactinomas = df[!, Symbol("Presence of hyperprolactinemia at diagnosis (0=no, 1=yes)")]
prolactinomas[ismissing.(prolactinomas)] .= 0
prolactinomas = Bool.(prolactinomas)

# prolactinomas = .~prolactinomas
# prolactinomas = df[prolactinomas.&&ividx, :]
prolactinomas = df[prolactinomas, :]
mean_prolactin = mean(prolactinomas[!, Symbol("If yes, please indicate prolactin levels (ng/ml)")])

low_prolactin = prolactinomas[!, Symbol("If yes, please indicate prolactin levels (ng/ml)")] .< mean_prolactin
low_prolactin = prolactinomas[low_prolactin, :]

high_prolactin = prolactinomas[!, Symbol("If yes, please indicate prolactin levels (ng/ml)")] .> mean_prolactin
high_prolactin = prolactinomas[high_prolactin, :]

# prolactinomas = high_prolactin
# prolactinomas = low_prolactin

# titlestring = "Prolactinomas, high prolactin"
# titlestring = "Prolactinomas, low prolactin"
titlestring = "Prolactinomas, all"
# titlestring = "Prolactinomas, no intervention"



growth = select(prolactinomas,
    [Symbol("PatientID"), Symbol("1st MRI date"),
        Symbol("2nd MRI date"),
        Symbol("3rd MRI date"),
        Symbol("4th MRI date"),
        Symbol("5th MRI date"),
        Symbol("6th MRI date"),
        Symbol("7th MRI date"),
        Symbol("8th MRI date"),
        Symbol("9th MRI date"),
        :X, :Y, :Z,
        :X_2, :Y_2, :Z_2,
        :X_3, :Y_3, :Z_3,
        :X_4, :Y_4, :Z_4,
        :X_5, :Y_5, :Z_5,
        :X_6, :Y_6, :Z_6,
        :X_7, :Y_7, :Z_7,
        :X_8, :Y_8, :Z_8,
        :X_9, :Y_9, :Z_9,
    ])


vol = zeros(9)
for i = 1:9
    if i == 1
        x = "X"
        y = "Y"
        z = "Z"
    else
        x = Symbol("X_" * string(i))
        y = Symbol("Y_" * string(i))
        z = Symbol("Z_" * string(i))

    end

    newvar = Symbol("Vol_" * string(i))
    growth[!, newvar] = 4 / 3 .* pi .* growth[!, x] ./ 2 .* growth[!, y] ./ 2 .* growth[!, z] ./ 2

end

results = DataFrame(PatientID=[], Time=[], Volume=[])
for (pt, i) in enumerate(growth[!, Symbol("PatientID")])
    begin
        timepts = Vector{Any}(undef, 9)
        timepts[1] = growth[pt, Symbol("1st MRI date")]
        timepts[2] = growth[pt, Symbol("2nd MRI date")]
        timepts[3] = growth[pt, Symbol("3rd MRI date")]
        timepts[4] = growth[pt, Symbol("4th MRI date")]
        timepts[5] = growth[pt, Symbol("5th MRI date")]
        timepts[6] = growth[pt, Symbol("6th MRI date")]
        timepts[7] = growth[pt, Symbol("7th MRI date")]
        timepts[8] = growth[pt, Symbol("8th MRI date")]
        timepts[9] = growth[pt, Symbol("9th MRI date")]
    end

    vols = Vector{Any}(undef, 9)
    vols[1] = growth[pt, :Vol_1]
    vols[2] = growth[pt, :Vol_2]
    vols[3] = growth[pt, :Vol_3]
    vols[4] = growth[pt, :Vol_4]
    vols[5] = growth[pt, :Vol_5]
    vols[6] = growth[pt, :Vol_6]
    vols[7] = growth[pt, :Vol_7]
    vols[8] = growth[pt, :Vol_8]
    vols[9] = growth[pt, :Vol_9]

    dfin = DataFrame(PatientID=fill(i, 9), Time=timepts, Volume=vols)

    sort!(dfin, :Time)

    t = zeros(Int, 9)
    for i = 1:length(timepts)
        wt = dfin[i, :Time] - dfin[1, :Time]
        println(wt)
        try
            t[i] = wt.value
        catch
            continue
        end
    end

    dfin.Time .= t
    sort!(dfin, :Time, rev=false)

    results = vcat(results, dfin)
end
dropmissing!(results)
filter!(:Volume => !=(0.0), results)

results.Time ./= 365
replace!(results.Volume, 0 => 10^-10)

CSV.write("parsed_data.csv", results)

# results[!, :PatientID] = convert(Vector{Int}, results[!, :PatientID])
# results[!, :PatientID] = convert(Vector{""}, results[!, :PatientID]),
results[!, :Time] = convert(Vector{Float64}, results[!, :Time])
results[!, :Volume] = convert(Vector{Float64}, results[!, :Volume])

#========================================================#
begin
    # Fit the model
    # m = fit(MixedModel, @formula(Volume ~ 1 + Time + (1 + Time | PatientID)), results)
    m = fit(MixedModel, @formula(Volume ~ Time + (Time | PatientID)), results)


    # Extract fixed-effect coefficients
    coefs = fixef(m)

    # Create a DataFrame for overall predictions
    n_times = 100
    new_time = range(minimum(results.Time), maximum(results.Time), length=n_times)
    pred_df = DataFrame(Time=new_time)

    # Calculate fixed-effects predictions manually
    pred_df[!, :Volume] = coefs[1] .+ coefs[2] .* pred_df.Time

    # Bootstrap to estimate uncertainty
    function bootstrap_prediction(data)
        m_boot = fit(MixedModel, @formula(Volume ~ 1 + Time + (1 + Time | PatientID)), data)
        coefs_boot = fixef(m_boot)
        return coefs_boot[1] .+ coefs_boot[2] .* pred_df.Time
    end

    bs = bootstrap(bootstrap_prediction, results, BasicSampling(1000))
    uncertainty = stderror(bs)

    # Plotting
    plt = plot(legend=:topleft)

    # Plot lines for individual patients
    for patient in unique(results.PatientID)
        patient_data = filter(row -> row.PatientID == patient, results)
        patient_data.Volume[patient_data.Volume .<= 0] .= 10^-10
        plot!(plt, patient_data.Time, patient_data.Volume, lw=2, color=:grey, alpha=0.5, legend=false)
    end

    # Add model prediction
    pred_df.Volume[pred_df.Volume .<= 0] .= 10^-10

    plot!(plt, new_time, pred_df.Volume, line=(2, :blue), label="Model Prediction", ribbon=uncertainty, legend=false,
        xlabel="Time [Year]", ylabel="Tumor volume [mm³]", dpi=300)
    # title!("Prolactinomas, no intervention")
    title!(titlestring)
    savefig("$titlestring log10.png")
    display(plt)
end
