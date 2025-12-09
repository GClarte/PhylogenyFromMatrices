## Phylogenies with Matricial Datasets

This is the implementation of the numerical methods described in my PhD to sample from the posterior associated with the model we propose to deal with Matricial datasets.

# "fonctions"

This folder contains all the functions needed to run the analysis. The code might be difficult to read.

# Reproducing the results of the papers

The R files at the root corresponds to the different experiments carried in the papers:
- Simulation.R corresponds to the first simulation without any misspecification
- Simulation_correlation.R to the simulation with correlated characters
- Simulation_iconicity1.R and iconicity2 to the two types of iconicity
- Simulation_tranformation.R to the unknown transformation
- Simulation_twotrees.R to the two trees dataset

- EuropeAsianLS.R corresponds to the presented result with Europe and Asian SL.
- The other files corresponds to other subparts of the dataset, most of which were not presented
- NJ.R computes raw representations of the initial dataset

On the Je languages:
- GenericTemplate_Je.R runs the simulations. The output of the study is stored in the AnalysisJeLanguages nexus files.
- The dataset is stored in je_languages_modified.csv

# Using the code

Everything is described in the GiveATry.R file, we did not include the data importation part, as this can be quite dependent on the datasets.
GiveAPlot.R plots all the parameters, and produces consensus trees, some parts might need to be adapted. The consensus tree must be annotated using for example treeannotator and then plotted with figtree.
