# CompariSAINT
Tool for formatting and comparing SAINTexpress outputs by fold change and t-test. 

This is a simple tool that accepts an unlimited number of CSV files holding SAINTexpress proteomics results with peptide counts. 

The script: 

1) Formats the output to have columns including: 

- Protein_ID
- Gene
- 2 x collapsed controls (Control_1 and Control_2)
- Sample spectral counts, marked by pool (biological replicate) and qe# (technical replicate)
- Average and sum spectral counts (AvgSpec and SpecSum)
- BFDR (Bayesian false discovery rate)

2) Concatenates the inputted SAINTexpress files to get one large dataframe allowing comparison of hits between baits

3) Finds the log2fc (log2 fold change) of the AvgSpec values for comparison, as well as the p-values, which are computed by scipy.stats.ttest_ind()

If you find this script useful, you are free to use it with attribution. 
