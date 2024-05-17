## Github Repository for
# Viral communities vary more across spatial scales than across annual sampling in soils recovering from decades long heating from an underground coal mine fire in Centralia, PA
## by Samuel Barnett and Ashley Shade
<i>This work is unpublished.</i>


### Data
Both the raw read data and metagenome assemblies for this study are available through NCBI under bioproject [PRJNA974462](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA974462/)

### To cite this work or code

TBD

### Abstract

Viruses are important components of the soil microbiome, influencing the population dynamics and function of their hosts. Yet, how soil viral populations are impacted by and recover from habitat disturbance, along with their changing influence on their microbial hosts, are still relatively unknown. With increasing attention to the long-term effects of climate and land use changes as well as other large scale anthropogenic environmental shifts, it is essential to understand how soil viruses recover from press disturbances. We used the unique conditions in Centralia, Pennsylvania, USA, to examine the dynamics of soil viral populations over seven consecutive years of recovery from decades-long soil heating. Underlying centralia is a coal mine fire that has been burning for over 60 years, significantly heating the soils above the fire front. As the fire moves along the coal seam, previously heated soils cool to ambient temperature allowing us to examine a gradient of heat disturbance intensity and recovery. We performed annual soil sampling and metagenomic sequencing to examining the viral populations across this gradient and as soils cooled over time. While we found significant changes in viral communities both over time and between fire affected and unaffected reference sites, we tended to see less resilience in the viral populations than we previous saw bacterial host community. Overall, dissimilarity in viral communities was greater when comparing across sites within a year relative to across years within a timepoint. We also measured significant changes in CRISPR investment as soils cooled, corresponding to shifts in viral diversity, together suggesting shifts in host selection for viral defense and predator-prey dynamics. Finally, we observed changes in viral-encoded auxiliary metabolic genes between fire affected and reference sites. Overall, these results indicate shifting host-viral interactions in soils undergoing recovery from long-term soil heating, with implications for soil community structure and function during and following recovery.

### Contents

Code is split up into two directories: [Sequence_processing](https://github.com/ShadeLab/Centralia_phages_Barnett/tree/main/Sequence_processing) and [Analysis](https://github.com/ShadeLab/Centralia_phages_Barnett/tree/main/Analysis).

#### Sequence processing
Code used for sequence processing including read QC, OTU clustering, taxonomy assignment, and tree building can be found under [Sequence_processing](https://github.com/ShadeLab/Centralia_phages_Barnett/tree/main/Sequence_processing). Scripts were run using SLURM on the MSU HPCC using slurm batch files with suffix .sb and are numbered by their order in the processing workflow. Outputs such as logs, warnings, or errors if any, are designated by the suffix .out and named in accordence with the library, run number, and slurm batch file. 

#### Analysis
Formal analysis can be found under [Analysis](https://github.com/ShadeLab/Centralia_phages_Barnett/tree/main/Analysis). All analysis was run with R and code was run in Rmarkdown. In the analysis directory you'll find the raw Rmarkdown files (.Rmd), a github friendly markdown rendering (.md) and the associated figure files from the rendering in separate sub-directories. The analysis was broken down into multiple chunks in separate Rmarkdown files:

### Funding
This work was supported by the U.S. National Science Foundation CAREER award 1749544. This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/).

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)
