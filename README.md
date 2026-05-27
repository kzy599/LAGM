# LAGM: Look-Ahead Genomic Mating

This repository contains the code accompanying the unpublished manuscript:
**"Look-ahead genomic mating optimizes the trade-off between genetic gain and inbreeding in family-based aquaculture breeding programs."**

## 📁 Project Structure

To reproduce the results reported in the study:

1. **Directory setup**

   Place all scripts in the **same directory**.

2. **Main script execution**

   Run the R simulation scripts:

   ```Shell
   Rscript Simulation_LAGM.r
   Rscript Simulation_LAGM_fixed.r
   ```

   These scripts generate all outputs used in the paper. The raw genomic data must be downloaded manually from NCBI. The preprocessing pipeline that converts the raw sequence data into the format consumed by the simulation lives in the [MNNDR](https://github.com/kzy599/MNNDR) repository.

3. **LAGM implementation**

   The LAGM method is implemented in the `lagmRcpp` subdirectory. Before running the simulations, install the package from GitHub:

   ```r
   remotes::install_github("kzy599/LAGM", subdir = "lagmRcpp")
   ```

4. **Dependencies**

   Before running any script, ensure that all required R libraries are available — they are loaded centrally via `utils.r` — and that **HiBLUP** (used for EBV calculation) is installed and configured correctly.

## Contact

If you encounter any issues or have questions about the code, please contact: `kangziyi1998@163.com`
