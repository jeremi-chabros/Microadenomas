using Plots
using CSV, DataFrames, XLSX, StatsBase, MixedModels, Dates, Bootstrap

df = DataFrame(XLSX.readtable("microadenoma_all.xlsx", 1))

intervention = df[!, Symbol("Any intervention during follow-up? (0=None, 1=Surgery, 2=Medical therapy, 3=Radiation therapy, 4=other)")]
# ividx = intervention .==0
# df = df[ividx, :]

pltxt = "All tumors, > median"

sizelist = [Symbol("1st Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("2nd Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("3rd Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("4th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("5th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("6th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("7th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("8th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("9th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("10th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("11th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("12th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("13th Recorded Tumor Size (A*B*C, in mm)"),
    Symbol("14th Recorded Tumor Size (A*B*C, in mm)")
]

growth = select(df, [
    :PatientID,
    Symbol("1st MRI date"),
    Symbol("2nd MRI date"),
    Symbol("3rd MRI date"),
    Symbol("4th MRI date"),
    Symbol("5th MRI date"),
    Symbol("6th MRI date"),
    Symbol("7th MRI date"),
    Symbol("8th MRI date"),
    Symbol("9th MRI date"),
    Symbol("10th MRI date"),
    Symbol("11th MRI date"),
    Symbol("12th MRI date"),
    Symbol("13th MRI date"),
    Symbol("14th MRI date"),
    sizelist...
])


dff = copy(growth)
for (seq, var) in enumerate(sizelist)
    df_sampl = copy(growth)
    replace!(df_sampl[!, var], missing => "")
    df_sampl[!, var] .= string.(df_sampl[!, var])

    function split_tmr(dimension)

        # Remove any non-numeric and non-"*" characters
        dimension = replace(dimension, r"[^\d*\.]" => "")

        # Split the string by "*"
        split_dim = split(dimension, "*")
        split_dim = filter(x -> x != "", split_dim)

        # Filter out empty strings
        split_dim = filter(x -> x != "", split_dim)

        # Return missing if no valid numbers are present
        if isempty(split_dim)
            return (x=missing, y=missing, z=missing)
        end

        # Convert the string values to floats
        float_dim = map(x -> parse(Float64, x), split_dim)


        # Compute x, y, z dimensions
        x = y = z = missing
        n = length(float_dim)

        if n == 1
            x = y = z = float_dim[1]
        elseif n == 2
            x, y = float_dim
            z = (x + y) / 2.0  # Average of x and y
        elseif n == 3
            x, y, z = float_dim
        else
            return (x=missing, y=missing, z=missing)
        end

        return (x=x, y=y, z=z)
    end

    vals = map(split_tmr, df_sampl[!, var])
    x = map(i -> vals[i].x, 1:length(vals))
    y = map(i -> vals[i].y, 1:length(vals))
    z = map(i -> vals[i].z, 1:length(vals))

    vol = 4pi / 3 .* x .* y .* z
    newvar = Symbol("Vol_$seq")
    dff[!, newvar] = vol
end


results = DataFrame(PatientID=[], Time=[], Volume=[])
for (pt, i) in enumerate(dff[!, Symbol("PatientID")])

    begin
        timepts = Vector{Any}(undef, 14)
        timepts[1] = dff[pt, Symbol("1st MRI date")]
        timepts[2] = dff[pt, Symbol("2nd MRI date")]
        timepts[3] = dff[pt, Symbol("3rd MRI date")]
        timepts[4] = dff[pt, Symbol("4th MRI date")]
        timepts[5] = dff[pt, Symbol("5th MRI date")]
        timepts[6] = dff[pt, Symbol("6th MRI date")]
        timepts[7] = dff[pt, Symbol("7th MRI date")]
        timepts[8] = dff[pt, Symbol("8th MRI date")]
        timepts[9] = dff[pt, Symbol("9th MRI date")]
        timepts[10] = dff[pt, Symbol("10th MRI date")]
        timepts[11] = dff[pt, Symbol("11th MRI date")]
        timepts[12] = dff[pt, Symbol("12th MRI date")]
        timepts[13] = dff[pt, Symbol("13th MRI date")]
        timepts[14] = dff[pt, Symbol("14th MRI date")]
    end

    vols = Vector{Any}(undef, 14)
    vols[1] = dff[pt, :Vol_1]
    vols[2] = dff[pt, :Vol_2]
    vols[3] = dff[pt, :Vol_3]
    vols[4] = dff[pt, :Vol_4]
    vols[5] = dff[pt, :Vol_5]
    vols[6] = dff[pt, :Vol_6]
    vols[7] = dff[pt, :Vol_7]
    vols[8] = dff[pt, :Vol_8]
    vols[9] = dff[pt, :Vol_9]
    vols[10] = dff[pt, :Vol_10]
    vols[11] = dff[pt, :Vol_11]
    vols[12] = dff[pt, :Vol_12]
    vols[13] = dff[pt, :Vol_13]
    vols[14] = dff[pt, :Vol_14]

    dfin = DataFrame(PatientID=fill(i, 14), Time=timepts, Volume=vols)

    sort!(dfin, :Time)

    t = zeros(Int, 14)
    for i = 1:length(timepts)
        wt = dfin[i, :Time] - dfin[1, :Time]
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
results.Time ./= 365

init_vol = results.Volume[results.Time.==0]
larger = results.Volume[results.Time.==0] .> median(init_vol)
smaller = results.Volume[results.Time.==0] .<= median(init_vol)

larger_id = results.PatientID[results.Time.==0][larger]
smaller_id = results.PatientID[results.Time.==0][smaller]

results = subset(results, :PatientID => x -> x .∈ Ref(larger_id))

filter!(:Volume => >=(0), results)
filter!(:Volume => <=(60_000), results)

results[!, :Time] = convert(Vector{Float64}, results[!, :Time])
results[!, :Volume] = convert(Vector{Float64}, results[!, :Volume])
# results.Volume[results.Volume.<=0] .= 10^-10

# CSV.write("parsed_data_all.csv", results)


#========================================================#
begin
    m = fit(MixedModel, @formula(Volume ~ Time + (Time | PatientID)), results)

    # fixed-effect coefficients
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

    bs_ci = confint(bs, PercentileConfInt(0.95))

    predictions = Float64[x[1] for x in bs_ci]
    lower_ci = Float64[x[2] for x in bs_ci]
    upper_ci = Float64[x[3] for x in bs_ci]

end

begin
    plt = plot(legend=:topleft)
    # lines for individual patients
    for patient in unique(results.PatientID)
        patient_data = filter(row -> row.PatientID == patient, results)
        # patient_data.Volume[patient_data.Volume.==10^-10] .= NaN
        plot!(plt, patient_data.Time, patient_data.Volume, lw=1, color=:grey, alpha=0.5, legend=false)
    end
    plot!(plt, new_time[2:end], predictions[2:end], line=(2, :blue), label="Model Prediction", legend=false,
        xlabel="Time [Year]", ylabel="Tumor volume [mm³]", dpi=300, ribbon=(upper_ci .- predictions, predictions .- lower_ci))

    title!(pltxt)
    savefig("plots/$pltxt.png")
    display(plt)
end


# # Log 10 plot
# begin
#     plt = plot(legend=:topleft)

#     # lines for individual patients
#     for patient in unique(results.PatientID)
#         patient_data = filter(row -> row.PatientID == patient, results)
#         y_data_log10 = log10.(patient_data.Volume)
#         y_data_log10[isinf.(y_data_log10)] .= NaN
#         plot!(plt, patient_data.Time, y_data_log10, lw=1, color=:grey, alpha=0.5, legend=false)
#     end

#     predictions_log10 = log10.(predictions[2:end])
#     predictions_log10[isinf.(predictions_log10)] .= NaN

#     upper_ci_log10 = log10.(upper_ci[2:end])
#     # lower_ci_log10 = log10.(lower_ci[2:end])
#     lower_ci_log10 = [lower_ci[i] .> 0 ? log10.(lower_ci[i]) : 0 for i in 2:length(lower_ci)]

#     plot!(plt, new_time[2:end], predictions_log10, line=(2, :blue), label="Model Prediction", legend=false,
#         xlabel="Time [Year]", ylabel="Log10(Tumor volume [mm³])", dpi=300, ribbon=(upper_ci_log10 .- predictions_log10, predictions_log10 .- lower_ci_log10))

#     plot!(plt, new_time[2:end], predictions_log10, line=(2, :blue), label="Model Prediction", legend=false,
#         xlabel="Time [Year]", ylabel="Log10(Tumor volume [mm³])", dpi=300, ribbon=log10.(uncertainty))

#     title!(pltxt)
#     savefig("plots/$pltxt log10.png")
#     display(plt)
# end
