# mvmr_simulations
Code to simulate data and run mvmr ancestry ajustment

Code is in 3 branches representing seperate simulation settings:

main- Main MVMR simualtions 
LD_mod_extra_run - Between sample heterogeneity simulations
realism_sims- MVMR sims using observed allele frequency patterns.

See specific README for details.


**MAIN**

**Input parameters**

Model – Define many of following parameters, refers to models_sims.R
Values
-A: B1=0, B2=0.8 	Simulates scenario where X1 does not effect Y, determining use of X2 in MVMR to adjust for false positive causal associations
          
-B: B1=0.4, B2=0.8	Both X1 and X2 have an effect on Y, X2 modifies magnitude of X-Y effect
          
-C: B1=0, B2=0		Neither X1 or X2 have an effect on Y
         
-D: B1=0.4, B2=0 	X1 has an effect on Y but X2 not associated
	
Pi – Varies X2->X1 effect, simulating correlation of exposure variable with ancestry 

SNPs/SNPsC   –  Varies number of SNPs generated for X1 and X2 respectively. Currently set as equal.   

Beta1/2 – Simulated causal effect of X1/ X2 on Y respectively. Set according to model (see above)

Reps – given number of times to loop through code. Results from each loop are saved 

Code structure
  
**  Modes_sims.R**
  	  Contains set up function to  define parameters in relation to the model passed in (A,B,C or D)
**  Functions_sims.R **
      Functions to perform most data generation/ analysis. Separated from main code for ease of interpretation
  
**  Run_sims.R**
  	  Runs simulations.
      Loops through each specified model in modes_sims.R and within each loop repeats for the given value of reps. Results for each loop and each model are saved together with indication of the run and model attached.
  
  
  
     
    General structure is as follows:
      For (model in A,B,C,D):
          	For rep in reps:
              Generate data (data_gen) – generates phenotype and genotype data based on the parameters defined in the model
              	                          GWAS on generated data (GWASres)
              Run two sample MR (univariate_MR) – uses TwoSampleMR package functions to perform separate two sample MR of X1->Y and X2-> Y
              Run multivariate MR (run_mvmr) – runs multivariate mr using the MVMR package
              

![image](https://github.com/user-attachments/assets/79d6b584-dfb0-4869-a1c0-d622c9357690)



  
  

