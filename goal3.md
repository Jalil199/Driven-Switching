Good. Now use the pilot runs to produce the first comparison and analysis layer.

Next task:

1. Run the small pilot parameter scan and collect all results.
2. Review which runs completed successfully and which ones failed or look suspicious.
3. Define or refine a clear and conservative switching criterion based on the pilot results.
4. Create a Jupyter notebook that automatically loads the pilot outputs and compares the runs. Use PyPlot for the plots with well chosen figures withe labels, units, etc... use this parameters for the plot 
rc("text", usetex=true)
rc("font", family="serif")
#rc("text.latex", preamble=raw"\usepackage{amsmath}")
rc("font", size=18)
rc("axes", labelsize=22, titlesize=22)
rc("xtick", labelsize=18)
rc("ytick", labelsize=18)
rc("legend", fontsize=16)
plt.rc("axes", linewidth=2)  # Set the linewidth of the plot axes.    
5. Include representative time traces, final magnetization comparisons, and a compact view of how the main observables change with parameter choice.
6. Clearly separate successful, failed, and suspicious runs in the analysis.
7. Comment all the function
Important constraints:

* Reuse the existing saved output structure.
* Keep the notebook readable and well commented.
* Focus on interpretability and correctness.
* Separete the jupyter notebook in different cells and add permitente coments
