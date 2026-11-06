**Phylogeography of the Endangered Franciscana Dolphin: Timing and Geological Setting of the Evolution of Populations**  

Code and input files for the biogeographic analysis (ancestral range estimation) of the franciscana dolphin (*Pontoporia blainvillei*).  
 
**Paper:** Nara L, Cremer M, Farro A, Colosio A, Barbosa L, Bertozzi C, Secchi E, Pagliani B, Costa-Urrutia P, Gariboldi M, Lazoski C, Cunha H (2022). Phylogeography of the Endangered Franciscana Dolphin: Timing and Geological Setting of the Evolution of Populations. *Journal of Mammalian Evolution*. [https://doi.org/10.1007/s10914-022-09607-7](https://doi.org/10.1007/s10914-022-09607-7 "https://doi.org/10.1007/s10914-022-09607-7")

**Overview**  
Ancestral geographic ranges were estimated with [BioGeoBEARS](http://phylo.wikidot.com/biogeobears "http://phylo.wikidot.com/biogeobears") on a time-stratified framework, comparing six models: DEC, DIVALIKE and BAYAREALIKE, each with and without the founder-event parameter *j*. Models were compared with likelihood ratio tests, AIC and AICc weights. 

**Repository contents**  
| | |  
|-|-|  
| **Path** | **Content** |   
| biogeography/biogeobears.R | Main analysis script (English, reference version) |   
| biogeography/pt/biogeobears_pt.R | Portuguese translation of the script, with explanatory comments |   
| biogeography/install_packages.R | Installs the R packages (run once) |   
| biogeography/input/ | Phylogeny, tip ranges, time periods and dispersal matrices (see input/README.md) |   
| biogeography/data/ | Accession number and localities of the samples from the study(see data/README.md) |   
| | | 

**Running the analysis**  
setwd("biogeography")  
 source("install_packages.R")   # once  
 source("biogeobears.R")  
   
On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/ "https://cran.r-project.org/bin/windows/Rtools/") before installing the packages.  
You might need to install some extra packages and dependencies.
If some of them were discontinued, you need to check the CRAN repository. 

**Data availability**  
The DNA sequences used in the study are available in GenBank. New haplotype sequences generated in this study are deposited under accession numbers OK188832, OK188833, OK188834 and OK188835. The full list of accessions is in data/accessions.tsv, and data/fetch_genbank.R downloads the sequences.  
The BioGeoBEARS input files (time-calibrated tree, tip ranges, time slices and dispersal matrices) were not archived at the time of publication. See biogeography/input/README.md.  

**Citation**  
If you use this code, please cite BioGeoBEARS:  
- Matzke NJ (2013). Probabilistic historical biogeography: new models for founder-event speciation, imperfect detection, and fossils allow improved accuracy and model-testing. *Frontiers of Biogeography* 5(4): 242-248.  
- Matzke NJ (2014). Model selection in historical biogeography reveals that founder-event speciation is a crucial process in island clades. *Systematic Biology* 63(6): 951-970.
  
**License and attribution**  
The analysis script is adapted from the example scripts distributed with BioGeoBEARS by N. J. Matzke and is released under GPL-3.0, as is BioGeoBEARS itself.  
