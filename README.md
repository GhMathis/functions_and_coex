########## TRIPLET COEXISTANCE, OSERMUNDIGEN 2022, LINK COEXISTANCE TO FUNCTIONS ##############
Data collection and code by Mathis Gheno on the project of Caroline Daniel. 
If you need informations ---> mathisgheno@gmail.com

files description:
- "code" : all the code made for the analyses.
	* "coex_data_oranisation" and "data_oranisation" : organize the data, not important
	* "enzymatic activity" : compute the contentration of p-nitrophenol release por B-gluosidase and Phosphatase
	* "growth_coef_computation" : compute 3 coef of growth rate for each plot
	* "link_function_and_coex" : main code and result in this file, plots are produce here. Link some functions
	to some coexistance values.

- "data" : all the data recolted, missing data : biomass_august_2022, enz_act_august_2022, cover_augsut
	* "biomass_OM_2022" : biomass, only spring 2022
	* "herbi_patho_OM_2022" : Herbivory and pathogen damage, spring and august
	* "par_and_height_OM_2022" : ligth interaction and canopy height, spring and august
	* "enzimatic_activity_OM_2022" : enzymatic activity, spring 
	* "biomcoexistence_all_2020_nitrogen_18sp2-3": coexistance data for all possible triplet with nitrogen
	* "final_triple_selection_corrected" : coexistance data for the 48 triplet, without nitrogen
	* "cover_OM_spring_2022" : cover data, spring

- "graph" : plot from the "link_function_and_coex" code file.

- "organized_data" : data in R.data format, cleaned and organezed for the analyses. This file contain the data
used in the code.

- "oster_data_2022" : data without modification (if needed). Can be removed.

- "workspaces" : the workspace, each one corresponding to a code file. Can be removed
