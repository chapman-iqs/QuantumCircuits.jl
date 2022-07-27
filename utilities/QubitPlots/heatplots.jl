function T2_ΓZ_plot(df::DataFrame)

	fids = df.sys_fids
	ΓZs = union(df.ΓZ);
	T2s = union(df.T2)

	fig, ax, hm = Makie.heatmap(ΓZs, T2s, reshape(fids, (4,3)), xlabel="ΓZ")
	ax.xlabel = "ΓZ (MHz)"
	ax.ylabel = "T2 (μs)"
	ax.title = "fidelities with target state"
	Makie.Colorbar(fig[:, end+1], hm)

	return fig
end
