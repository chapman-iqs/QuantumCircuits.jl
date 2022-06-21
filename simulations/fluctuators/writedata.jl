function exportpath(path_to_data, foldername)
	ep = string(path_to_data, foldername, "/", Dates.now(), "/")
	replace(ep, ":" => "-")
end

function writedata(pars::OrderedDict, dataframes::Vector{DataFrame}, datanames::Vector{String},
					path_to_data = datapath, foldername = defaultfolder)

	ep = exportpath(path_to_data, foldername)
	mkpath(ep)
	println(string("Writing to path ", ep))

	CSV.write(string(ep, "parameters.csv"), pars)

	for (data, name) in zip(dataframes, datanames)
		CSV.write(string(ep, name, ".csv"), data)
	end
end

function writedata(pars::OrderedDict, dataframe::DataFrame, dataname::String,
					path_to_data = datapath, foldername = defaultfolder)

	ep = exportpath(path_to_data, foldername)
	mkpath(ep)
	println(string("Writing to path ", ep))

	CSV.write(string(ep, "parameters.csv"), pars)
	CSV.write(string(ep, dataname, ".csv"), dataframe)
end
