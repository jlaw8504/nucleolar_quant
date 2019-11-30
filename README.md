# nucleolar_quant
MATLAB functions that measure fluorescently labeled nucleolar signals. Jupyter notebooks summarizing nucleolar signal data as seaborn violin plots are also provided in the jupyter_notebooks folder.
## MATLAB (.m files) Requirements
- MATLAB
- MATLAB Image Processing Toolbox
- MATLAB Bioformats Plugin - available at https://www.openmicroscopy.org/bio-formats/downloads/
## Python/Jupyter Notebook requirements
- JSON
- matplotlib
- seaborn
- pandas
# MATLAB Code Organization
* rdna_corr_summary.m
   * rdna_corr.m
* cdc14_cbf5_compaction_summary.m
  * quant_int.m
* tetra_dimer_summary.m
  * quant_int.m
* mutant_summary.m
  * mutant_analysis.m
   * quant_int
* pNOY130_area_summary.m
* pNOY_plasmid_summary.m
  * dir_rhos
   * parse_nd2
    * MATLAB Bioformats Plugin
   * corr_cdc14_cbf5
# Jupyter Notebooks
The Jupyter Notebooks parse the JSON files outputted by MATLAB summary programs to generate seaborn violin plots.
