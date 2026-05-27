# LAGM: Look-ahead genomic mating


This repository contains the code used in the unpublished manuscript:  
**"Look-ahead genomic mating optimizes the trade-off between genetic gain and inbreeding in family-based aquaculture breeding programs"**

## 📁 Project Structure

To reproduce the results from the study:

1. **Directory Setup**  

   Ensure all scripts are placed in the **same directory**.

3. **Main Script Execution**  

   Run the R script:

   ```Shell
   Rscript Simulation_LAGM.r
   Rscript Simulation_LAGM_fixed.r
   ```
   This will generate all outputs used in the paper. The raw genomic data must be manually downloaded from NCBI. The preprocessing code for preparing the raw sequence data before inputting it into the simulation can be found in the MNNDR repository.

4. **LAGM Implementation**

   see in lagmRcpp, before runing, install the packages with remotes::install_github("kzy599/LAGM", subdir = "lagmRcpp")

   
5. **Dependencies**

   Before running any scripts, make sure all required libraries of  R `utils.r` and software for EBV calculationg (HiBlup) are installed correctly.

## Contact

If you encounter any issues or have questions about the code, feel free to contact me at:`kangziyi1998@163.com`
