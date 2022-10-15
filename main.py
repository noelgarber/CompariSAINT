import numpy as np
import pandas as pd
import math
import warnings
import scipy.stats

alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

#Create a dictionary of baits to compare

bait_df_dict = {}

no_more_baits = False

while not no_more_baits: 
	bait_filename = input("Enter a bait filename (.csv) for comparison, or hit enter when done:  ")
	if bait_filename != "": 
		bait_df = pd.read_csv(bait_filename)
		bait_name = bait_df.at[1, "Bait"] #Takes the first bait name, which should be repeated/the same for the whole column
		bait_df_dict[bait_name] = bait_df
	else: 
		no_more_baits = True

#----------------------------------------------------------------------------------------------------------------------------

#Format each bait dataframe

keep_columns = ["Prey", "PreyGene", "ctrlCounts", "Spec", "AvgSpec", "SpecSum", "BFDR"]

reps_cols_dict = {} #Empty dict that will contain lists of tech reps lists for each bio rep for each bait
bait_cols_dict = {} #Empty dict that will contain the column names for bait data

def maxima(values_string, separator = "|", entries = 2): 
	values = values_string.split(separator)
	values = [eval(i) for i in values]
	values.sort(reverse = True) #sorts descending
	output_list = []
	for i in np.arange(entries):
		output_list.append(values[i])
	return output_list

def split_data(values_string, separator = "|"): 
	values = values_string.split(separator)
	values = [eval(i) for i in values]
	return values

print("Formatting dataframes...")

for bait_name, bait_df in bait_df_dict.items(): 
	bait_df = bait_df[keep_columns] #Keeps only the desired columns

	bait_df["ctrlCounts"] = bait_df["ctrlCounts"].apply(maxima) #Extract the top controls
	bait_df["Spec"] = bait_df["Spec"].apply(split_data)

	#Make column names to hold spectral counts
	bio_reps = int(input("\tFor " + bait_name + ", how many biological replicates? Enter an integer:  "))
	tech_reps = int(input("\tFor " + bait_name + ", how many technical replicates per biological replicate? Enter an integer:  "))
	
	reps_lists = []
	spec_cols = []

	for i in np.arange(1, bio_reps + 1): 
		bio_reps_list = []
		for j in np.arange(1, tech_reps + 1): 
			col = bait_name + "_Pool_" + alphabet[i - 1] + "_qe" + str(j)
			spec_cols.append(col)
			bio_reps_list.append(col)
		reps_lists.append(bio_reps_list)

	#Split the top controls and spectra into separate columns of individual values
	controls_split_df = pd.DataFrame(bait_df["ctrlCounts"].tolist(), columns = ["Control_1", "Control_2"]) #Make dataframe with split controol values
	spec_split_df = pd.DataFrame(bait_df["Spec"].tolist(), columns = spec_cols) #Make dataframe with split spectral counts
	bait_df = pd.concat([bait_df, controls_split_df, spec_split_df], axis=1) #Concatenate split values back into original dataframe
	bait_df = bait_df.drop(["ctrlCounts", "Spec"], axis=1) #Drop the original unsplit values, as they are no longer needed

	#Rename specific columns to be more descriptive
	col_rename_dict = {
		"Prey": "Protein_ID", 
		"PreyGene": "Gene", 
		"AvgSpec": bait_name + "_AvgSpec", 
		"SpecSum": bait_name + "_SpecSum", 
		"BFDR": bait_name + "_BFDR"
	}

	bait_df.rename(columns = col_rename_dict, inplace = True)

	#Set index as protein ID
	bait_df.set_index("Protein_ID", drop = False, inplace = True)

	#Reassign back to dict
	bait_df_dict[bait_name] = bait_df

	#Make list of bait-specific cols
	cols = spec_cols
	cols.extend([bait_name + "_AvgSpec", bait_name + "_SpecSum", bait_name + "_BFDR"])
	bait_cols_dict[bait_name] = cols #Assign to dict
	reps_cols_dict[bait_name] = reps_lists

	print("\tDone!")

#----------------------------------------------------------------------------------------------------------------------------

#Concatenate dataframes

print("Concatenating dataframes...")

#Create non-redundant list of hits present in at least one of the bait dataframes
hits_list = []

for bait_name, bait_df in bait_df_dict.items(): 
	#Add hits to hits_list
	hits = bait_df["Protein_ID"].tolist()
	hits_list.extend(hits)

hits_list = list(dict.fromkeys(hits_list)) #Removes duplicates

all_baits_df = pd.DataFrame(index = hits_list) #No columns declared yet, these will be declared dynamically
all_baits_df.index.name = "Protein_ID"

#Construct concatenated dataframe

for protein_id in all_baits_df.index.values.tolist(): 
	#Iterate through bait dataframes
	for bait_name, bait_df in bait_df_dict.items(): 
		bait_specific_cols = bait_cols_dict.get(bait_name)
		cols = ["Gene", "Control_1", "Control_2"]
		cols.extend(bait_specific_cols)
		lookup_results = {}
		for col in cols: 
			try: 
				all_baits_df.at[protein_id, col] = bait_df.at[protein_id, col]
			except KeyError: 
				if col != "Gene" and "BFDR" not in col: 
					all_baits_df.at[protein_id, col] = 0
				elif "BFDR" in col: 
					all_baits_df.at[protein_id, col] = 1

print("\tDone!")
print("---")

#----------------------------------------------------------------------------------------------------------------------------

#Compare AvgSpec values with log2fc and t-tests

bait_list = list(bait_df_dict.keys())

#Permute list of baits to get pairs for comparison
bait_pairs = []
permute_pairs = input("For bait-bait comparisons, permute pairs? If no, manual input will be prompted. (Y/N)  ")

if permute_pairs == "Y": 
	for bait1 in bait_list: 
		for bait2 in bait_list: 
			if bait1 != bait2: 
				pair = (bait1, bait2)
				bait_pairs.append(pair)
else: 
	no_more_pairs = False
	print("bait_list =", bait_list)
	while not no_more_pairs: 
		bait1 = input("Enter bait #1 of a pair to compare, or hit enter if done:  ")
		bait2 = input("Enter bait #2 of the pair (or hit enter to abort to last pair inputted):  ")
		if bait1 != "" and bait2 != "": 
			bait_pairs.append((bait1, bait2))
		else: 
			no_more_pairs = True
		print("---")

#Define a log2fc function that takes two values, as well as a value to add to avoid mathematical errors (since log(0) is invalid)
def log2fc(value1, value2, add = 0.1): 
	log_val1 = math.log2(value1 + add)
	log_val2 = math.log2(value2 + add)
	log2fc = log_val2 - log_val1
	return log2fc

#Define a function to return lists of tech-rep-mean values for biological replicates
def bio_reps_means(df_index, reference_df, col_labels_tuple): 
	list_out = []
	for bio_rep in col_labels_tuple: 
		tech_rep_values = []
		for tech_rep in bio_rep: 
			tech_rep_value = reference_df.at[df_index, tech_rep]
			tech_rep_values.append(tech_rep_value)
		bio_rep_mean = sum(tech_rep_values) / len(tech_rep_values)
		list_out.append(bio_rep_mean)
	return list_out

#Evaluate log2fc and t-test results
for protein_id in all_baits_df.index.values.tolist(): 
	for bait_pair in bait_pairs: 
		#Compute and assign log2fc
		bait1, bait2 = bait_pair
		avgspec1, avgspec2 = all_baits_df.at[protein_id, bait1 + "_AvgSpec"], all_baits_df.at[protein_id, bait2 + "_AvgSpec"]
		fc = log2fc(avgspec1, avgspec2)
		all_baits_df.at[protein_id, bait1 + "_" + bait2 + "_log2fc"] = fc

		#Compute and assign p-value from t test of biological replicates
		bait1_reps_cols, bait2_reps_cols = reps_cols_dict.get(bait1), reps_cols_dict.get(bait2)

		bait1_bio_values = bio_reps_means(df_index = protein_id, reference_df = all_baits_df, col_labels_tuple = bait1_reps_cols)
		bait2_bio_values = bio_reps_means(df_index = protein_id, reference_df = all_baits_df, col_labels_tuple = bait2_reps_cols)

		if len(bait1_bio_values) < 3: 
			warnings.warn("\tThere are less than 3 biological replicates for " + bait1 + "; t-test results using it will have low confidence.")
		elif len(bait2_bio_values) < 3: 
			warnings.warn("\tThere are less than 3 biological replicates for " + bait2 + "; t-test results using it will have low confidence.")

		t_statistic, pvalue = scipy.stats.ttest_ind(bait1_bio_values, bait2_bio_values, equal_var = False)

		all_baits_df.at[protein_id, bait1 + "_" + bait2 + "_welchs_pvalue"] = pvalue

#Sort by best BFDR

bfdr_cols = []
for bait in bait_list: 
	bfdr_cols.append(bait + "_BFDR")

def best_bfdr_getter(index, dataframe): 
	bfdr_list = []
	for bfdr_col in bfdr_cols: 
		bfdr_list.append(dataframe.at[index, bfdr_col])
	min_bfdr = min(bfdr_list)
	return min_bfdr

for protein_id in all_baits_df.index.values.tolist(): 
	best_bfdr = best_bfdr_getter(index = protein_id, dataframe = all_baits_df)
	all_baits_df.at[protein_id, "Best_BFDR"] = best_bfdr

all_baits_df.sort_values("Best_BFDR", ascending = True, inplace = True)
all_baits_df.drop("Best_BFDR", axis = 1, inplace = True)

#Move birA and SAV_STRAV to the top

for protein_id in all_baits_df.index.values.tolist(): 
	gene = all_baits_df.at[protein_id, "Gene"]
	if gene == "birA": 
		all_baits_df.at[protein_id, "Sorter"] = 2
	elif gene == "SAV_STRAV": 
		all_baits_df.at[protein_id, "Sorter"] = 1
		all_baits_df.at[protein_id, "Gene"] = "Streptavidin"

all_baits_df.sort_values("Sorter", ascending = False, inplace = True)
all_baits_df.drop("Sorter", axis = 1, inplace = True)

print("Final dataframe:")
print(all_baits_df)

#Save output
all_baits_df.to_csv("comparisaint_output.csv")