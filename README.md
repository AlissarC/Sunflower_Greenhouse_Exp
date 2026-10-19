# Drought × Nitrogen Sunflower Greenhouse Experiment Data Repository

This repository contains data and R scripts from a drought-by-nitrogen sunflower greenhouse experiment conducted in the [Smith Ecophys Lab](http://www.smithecophyslab.com/) at Texas Tech University.

The data and analyses in this repository are associated with the manuscript:

**"Water availability alters water and nutrient acquisition costs: Implications for plant resource economics"**  
Cheaib & Smith (2026), submitted to *Annals of Botany*.

## Repository contents

This repository contains the data used for the results presented in the manuscript, associated metadata describing the variables included in the data files, and the R scripts used to process and analyze the data and generate the figures and tables.

### `input/`

The `input` folder contains the data used in the analyses, including:

- Transpiration data
- Biomass data
- Leaf trait data, including gas exchange (A-Ci curves), leaf nitrogen, and isotope measurements
- Plant morphology data
- Leaf and leaf-disk scan images
- Metadata files describing the variables, units, and other relevant information associated with the datasets 

### `scripts/`

The `scripts` folder contains the R scripts used for data processing, statistical analyses, and generation of figures and tables presented in the manuscript.

Main scripts include:

- `scripts/Sunflower_processing.R` — data processing and calculation of traits and variables
- `scripts/Stat_sunflower_t4_final.R` — statistical analyses
- `scripts/Fig_pannels.R` — generation of figure panels
- `scripts/tables.R` — generation of tables

Figures were generated using these R scripts, except for **Figure 1b and Figure 5**, which were created using BioRender.

### `scripts/functions/`

The `scripts/functions/` folder contains R functions called by the main scripts, including functions used to calculate leaf-level carbon costs.

### `output/`

The `output` folder contains final files used in the analysis.

### `Fig/`

The `Fig` folder contains final figures presented in the manuscript.

### `tables/`

The `tables` folder contains the outputs of the lmms applied and presented in the manuscript.


## Citation

If you use these data or scripts, please cite:

Cheaib & Smith (2026). *Water availability alters water and nutrient acquisition costs: Implications for plant resource economics*. Submitted to *Annals of Botany*.

