##### BioGeoBEARS: ancestral range estimation (DEC, DIVALIKE, BAYAREALIKE, each with and without +J) #####
# Adapted from the BioGeoBEARS example scripts by N. J. Matzke
# (http://phylo.wikidot.com/biogeobears; Matzke 2013, Frontiers of Biogeography 5:242-248; Matzke 2014, Systematic Biology 63:951-970)
# Adapted for: Phylogeography of the Endangered Franciscana Dolphin: Timing and Geological Setting
# of the Evolution of Populations. J Mammal Evol (2022). https://doi.org/10.1007/s10914-022-09607-7
# Changes: study-specific input files and settings (see input/README.md)
# Packages: run install_packages.R once. Run this script from the biogeography/ folder. License: GPL-3

dir.create("results", showWarnings=FALSE)

##### Script BioGeoBEARS #####

# To install BioGeoBEARS you need to install the Rtool

# Install:


# Open:

library(Rcpp)
library(RcppArmadillo)
library(gdata)
library(gtools)
library(xtable)
library(plotrix)
library(vegan)
library(FD)
library(SparseM)
library(ape)
library(rexpokit)
library(cladoRcpp)
library(GenSA)
library(hms)
library(phytools)
library(expm)
library(fdrtool)
library(statmod)
library(spam)
library(MultinomialCI)
library(optimx)
library(snow)
library(phylobase)
library(devtools)

#BioGeoBEARS is installed using devtools and directly from the github:


# The package is developed to ML and Bayesian approaches
# (a) ancestral geographic ranges, and
# (b) models for the evolution of geographic range across a phylogeny.

# Open:

library(BioGeoBEARS)

# You will need to set your working directory:

# Run this script from the biogeography/ folder

# Extra files are stored in the "extdata" directory inside the installed package.
# BioGeoBEARS contains various example files in its extdata directory.
# To find the "exdata" directory:

extdata_dir = normalizePath(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# Indicate where the phylogeny and geographic distribution files are:

trfn = "input/tree.nwk"
geogfn = "input/geography.data"

# Phylogeny file #
#  1.  Must be binary/bifurcating: no polytomies
#  2.  No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
#  3.  Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
#  4.  You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the top of the tree as fossils.
#  5.  The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of
#      millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#      e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#      your branchlengths to get them into reasonable units.
#  6.  Don't use spaces in the species names, use E.G. "_"

# Check the raw Newick file:

moref(trfn)

# Check the phylogeny

pdffn = "results/tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr)
title("TITLE")
axisPhylo(side = 1, root.time = NULL, backward = TRUE)

dev.off()

# Geographic file #
#  1.   It is a PHYLIP-formatted file.
#       This means that in the first line:
#    -  1st number is the number of lines (species/demes)
#    -  2nd number is the number o columns (number of areas)
#    -  after a tab, put the areas in parentheses, with spaces: (A B C D)
#
#  1.5. Example first line:
#    10    4    (A B C D)
#
#  2.   The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
#
#  2.5a. This means a TAB between the species name and the area 0/1s
#  2.5b. This also means NO SPACE AND NO TAB between the area 0/1s
#  3.    All names in the geographic files must match the names in the phylogeny file
#  4.    Don't use spaces in the species names, use E.G. "_"
#  5.    Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, i.e. genetically isolated populations.
#    These may or may not be identical with species. You would NOT want to just use specimens, as each specimen
#    automatically can only live in 1 area, which will typically favor DEC+J. This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. In summary, you should collapse multiple specimens into species/lineages if
#    data indicates they are the same genetic population.

# Check the raw geography file:
moref(geogfn)

# Check the geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger
# than the number of areas you set up, but it can be smaller.

max_range_size = 12

# KEY HINT: The number of states (= number of different possible geographic ranges)
# depends on (a) the number of areas and (b) max_range_size. If you have more than
# about 500-600 states, the calculations will get REALLY slow, since the program has
# to exponentiate a matrix of e.g. 600x600.  Often the computer will just sit there and
# crunch, and never get through the calculation of the first likelihood.
# To check the number of states for a given number of ranges, try:
#numstates_from_numareas(numareas=12, maxareas=9, include_null_range=TRUE)

# Large numbers of areas have problems ...unless you limit the max_range_size:
#numstates_from_numareas(numareas=12, maxareas=4, include_null_range=TRUE)

#################################################################################
# DEC and DEC+J #
#################################################################################
# The BioGeoBEARS "DEC" model is identical with the Lagrange DEC model, and should return
# identical ML estimates of parameters, and the same log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly different, since BioGeoBEARS is
# reporting the ancestral state probabilities under the global ML model, and Lagrange is
# reporting ancestral state probabilities after re-optimizing the likelihood after fixing
# the state at each node. These will be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.

# Intitialize a default model (DEC model):
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file:
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file:
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size:
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min to treat tip as a direct ancestor (no speciation event):
# BioGeoBEARS_run_object$min_branchlength = 0.000001

# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = TRUE

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree...", below.
# Uncomment files you wish to use in time-stratified analyses:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing, if desired:
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object.
# Required for these to work!(It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states:
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # obter os estados ancestrais da corrida optim

# Set up DEC model(nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just load the saved result.

runslow = TRUE
resfn = "results/DEC_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

#################################################################################################################################
# Run DEC+J
#################################################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis):
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table::
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states:
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results/DEC+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resDECj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
    }


# PDF plots
pdffn = "results/DEC+J_v1.pdf"
pdf(pdffn, width=6, height=6)

###############################################################################################
# Plot ancestral states - DEC #
###############################################################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

###############################################################################################
# PPlot ancestral states - DECJ #
###############################################################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF

##############################################################################################################
# DIVALIKE and DIVALIKE+J #
##############################################################################################################
# The BioGeoBEARS "DIVALIKE" model is not identical with Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's processes the way DEC does, but only allowing the
# processes DIVA allows (widespread vicariance: yes; subset sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.
##############################################################################################################

##############################################################################################################
# Run DIVALIKE #
##############################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Sem a simpatrica com subset
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "results/DIVALIKE_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resDIVALIKE = res
    } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
    }

#######################################################
# Run DIVALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Obtem os valores de parametros ML do modelo nested de 2-parametros
# (isso garante que o modelo de 3-parametros seja sempre ao menos bom)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results/DIVALIKE+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")

    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resDIVALIKEj = res
    } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
    }

pdffn = "results/DIVALIKE_vs_DIVALIKE+J_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

###############################################################################################
# BAYAREALIKE and BAYAREALIKE+J #
###############################################################################################
# As with DIVA, the BioGeoBEARS BayArea-like model is not identical with the full Bayesian model
# implemented in the "BayArea" program of Landis et al. (2013).
#
# Instead, this is a simplified likelihood interpretation of the model.  Basically, in BayArea and
# BioGeoBEARS-BAYAREALIKE, "d" and "e" work like they do in the DEC model of Lagrange (and BioGeoBEARS),
# and then BayArea's cladogenesis assumption (which is that nothing in particular happens at cladogenesis)
# is replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is to test the importance of the cladogenesis model on
# particular datasets. Does it help or hurt the data likelihood if there is no special cladogenesis process?
#
# BAYAREALIKE is a likelihood interpretation of BayArea, and it is "like BayArea" -- similar to, but not
# identical to, Bayesian BayArea.
###############################################################################################

###############################################################################################
# Run BAYAREALIKE
###############################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# Sem subset de simpatria
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "results/BAYAREALIKE_v1.Rdata"
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)
    resBAYAREALIKE = res
    } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
    }

###################################################################################################################
# Run BAYAREALIKE+J
###################################################################################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Set up a time-stratified analysis:
BioGeoBEARS_run_object$timesfn = "input/timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "input/manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "input/areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "input/areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "input/distances_matrix.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$num_cores_to_use = 1

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# BAYAREALIKE+J seems to crash on some computers, usually Windows
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max"
# slightly for each parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "results/BAYAREALIKE+J_v1.Rdata"
runslow = TRUE
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res

    save(res, file=resfn)

    resBAYAREALIKEj = res
    } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
    }

pdffn = "results/BAYAREALIKE_vs_BAYAREALIKE+J_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
analysis_titletxt ="TITLE"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()

#######################################################################
# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################################
# Statistics - DEC vs. DEC+J
#######################################################################
# We have to extract the log-likelihood differently, depending on the version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models confer the same likelihood on the data.

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
# We have to extract the log-likelihood differently, depending on the version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################################
# ASSEMBLE RESULTS TABLES:  DEC, DEC+J, DIVALIKE, DIVALIKE+J, ...     #
#######################################################################
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# To see the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for e.g. convergence issues
#######################################################

# Loads to "restable"
save(restable, file="results/restable_v1.Rdata")
load(file="results/restable_v1.Rdata")

# Loads to "teststable"
save(teststable, file="results/teststable_v1.Rdata")
load(file="results/teststable_v1.Rdata")

# Save to text file
write.table(restable, file="results/restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="results/teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs – factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

# Save in text file
write.table(restable_AIC_rellike, file="results/restable_AIC_rellike.txt", quote=FALSE, sep="\t")
write.table(restable_AICc_rellike, file="results/restable_AICc_rellike.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike), file="results/restable_AIC_rellike_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc_rellike), file="results/restable_AICc_rellike_formatted.txt", quote=FALSE, sep="\t")

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
