# Generation-Adaptive Selection for Balancing Genetic Gain and Diversity in Family-Based Aquaculture Breeding Programs

This repository contains the code used in the unpublished manuscript:  
**"Generation-Adaptive Selection for the Balance Between Genetic Gain and Diversity in Family-Based Aquaculture Breeding Programs."**

## 📁 Project Structure

To reproduce the results from the study:

1. **Directory Setup**  

   Ensure all R and Python scripts are placed in the **same directory**.

3. **Main Script Execution**  

   Run the R script:

   ```Shell
   Rscript main.r
   ```
   This will generate all the outputs used in the paper.This will generate all outputs used in the paper. The raw genomic data must be manually downloaded from NCBI. The preprocessing code for preparing the raw sequence data before inputting it into the simulation can be found in the MNNDR repository.

4. **GAS Implementation**

   The core implementation of the Generation-Adaptive Selection (GAS) algorithm is located in the Python script `optMatingP.py`.
   
5. **Dependencies**

   Before running any scripts, make sure all required libraries of both R `utils.r` and Python `optMatingP.py` are installed correctly.

## Contact

If you encounter any issues or have questions about the code, feel free to contact me at:`kangziyi1998@163.com`
