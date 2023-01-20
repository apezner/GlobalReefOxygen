# *Increasing hypoxia on global coral reefs under warming*

GitHub repository containing code accompanying the global coral reef oxygen manuscript (Pezner et al., 2023; DOI: [TBD])

**Data Repository DOI:** [TBD]

**Authors:** Ariel K. Pezner, Travis A. Courtney, Hannah C. Barkley, Wen-Chen Chou, Hui-Chuan Chu, Samantha M. Clements, Tyler Cyronak, Michael D. DeGrandpre, Samuel A. H. Kekuewa, David I. Kline, Yi-Bei Liang, Todd R. Martz, Satoshi Mitarai, Heather N. Page, Max S. Rintoul, Jennifer E. Smith, Keryea Soong, Yuichiro Takeshita, Martin Tresguerres, Yi Wei, Kimberly K. Yates, and Andreas J. Andersson

**Abstract:** Ocean deoxygenation is predicted to threaten marine ecosystems globally. However, current and future oxygen concentrations and the occurrence of hypoxic events on coral reefs remain underexplored. Here, using autonomous sensor data to explore oxygen variability and hypoxia exposure at 32 representative reef sites, we reveal that hypoxia is already pervasive on many reefs. 84% of reefs experienced weak to moderate (≤153 to ≤92 µmol O<sub>2</sub> kg<sup>-1</sup>) hypoxia and 13% experienced severe (≤61 µmol O<sub>2</sub> kg<sup>-1</sup>) hypoxia. Under different climate change scenarios based on 4 Shared Socioeconomic Pathways (SSPs), we show that projected ocean warming and deoxygenation will increase the duration, intensity, and severity of hypoxia, with more than 94% and 31% of reefs experiencing weak to moderate and severe hypoxia, respectively, by 2100 under SSP5-8.5. This projected oxygen loss could have negative consequences for coral reef taxa due to the key role of oxygen in organism functioning and fitness.


**Citation:** [TBD]

---

### Repository contains the following:

1. [code](https://github.com/apezner/GlobalReefOxygen/tree/master/code)
  * ***PeznerMS.R*** - R Code used to analyze the dissolved oxygen (DO) dataset under present-day conditions and calculate projected changes in DO under 5 warming scenarios for each coral reef site (creates Figures 2-3, Extended Data Figures 1-2 and 6-8, and statistics); requires ***calcDOatsat.R***.
  * ***CMIP6_loc.R*** - R Code used to extract temperature projection data from the Coupled Model Intercomparison Project 6 (CMIP6) ensemble member Community Earth System Model Whole Atmosphere Community Climate Model (CESM2-WACCM) model for each of the 12 coral reef locations (creates Extended Data Figure 3). CESM2-WACCM model data used may be downloaded [here](https://doi.org/10.22033/ESGF/CMIP6.10028) and [here](https://doi.org/10.22033/ESGF/CMIP6.10101).
  * ***calcDOatsat.R*** - R Function used to calculate dissolved oxygen solubility in seawater (required for ***PeznerMS.R***)
  * ***plotmap_multiplot_figure.R*** - R Code used to map coral reef locations, calculate DO climatologies for all locations, and calculate distributions of DO for all locations (creates Figure 1).
  * ***boxmodel.m*** - MATLAB code used to create and run a simple coral reef box model to assess calculation approach used in ***PeznerMS.R*** (creates Extended Data Figure 5); requires ***oxy_sol.m*** and ***reefo2dif.m*** functions.
  * ***oxy_sol.m*** - MATLAB function used to calculate dissolved oxygen solubility in seawater (required for ***boxmodel.m***)
  * ***reefo2dif.m*** - MATLAB function with differential equations required by box model (required for ***boxmodel.m***)
