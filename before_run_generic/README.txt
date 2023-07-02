This folder contains all the files necessary to replicate "Explaining Intergenerational Mobility: The Role of Fertility and Family Transfers" by Daruich and Kozlowski.
The codes make use of the Compecon library available at http://www4.ncsu.edu/~pfackler/compecon/toolbox.html. Please make sure that the Compecon library is properly compiled following the instructions on the website. We also use a cluster to solve the simulations of the model.

- The initial folder contains all the master files to do all the exercises:
	- start_calibrate.m solves the model, simulates and computes the moments of interest
	- run_compare_ss.m computes the cross-state validation exercises of Table 3
	- run_VLEIC_age.m computes the variance of future earnings decomposition used in Figure 4. Note that to run this code we used a computer with 250Gb of RAM. Figure_VLEIC.m produces such figure.
	- run_model_options.m computes the moments needed for Tables 6, 7 and 8. Note that to run this code we used a computer with 250Gb of RAM.
		- Fig_Tables_Models_Options_Slides.m summarizes the results used in Tables 6 and 7.  
		- Fig_Tables_Models_Options_Policy.m summarizes the results used in Table 8.  
		
- Folder \run_generic\ contains all subfiles used to solve the model, simulate and calculate moments.
- Folder \results\ stores all the results.
	- In here we already provide the results we obtained using these codes.