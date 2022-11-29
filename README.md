# climatic_debt

New project on climatic debt in planktonic foraminifera.

## File overview:

**config_file.R**
- general configurations for visualisations

**create_depth_tbl.R**
- creates a table with preferred water depth of individual foraminifera species

**clean_data.R**
- cleans species level data based on taxonomy and occurrences  
  
**extract_temperature.R**
- adds the ambient temperature to each species occurrence at the surface and the preferred depth based on gcm models  
  
**calculate_climatic_debt.R**
- calculates the mismatch between preferred temperature of assemblages and actual temperature on the surface  
  
**analyze_climatic_debt.R**
- calculates general trends of climatic debt per latitudinal zone, as well as range debt  
  
**model_climatic_debt.R**
- calculates climatic debt through time and links it to lags in temperature  
  
**robustness_depth.R**
- robustness test calculating climatic debt based on temperature in the preferred depth of focal species instead of surface temperature  
  
- **robustness_abundance.R**
- robustness test calculating climatic debt based on presence-absence of species instead of relative abundance

- temperature proxy data was taken from [Friedrich & Timmermann 2020](https://doi.org/10.1016/j.epsl.2019.115911)