# PhyloNets
Create Phylogenetic Networks from Clusters

This repository contains the source files for the results presented in the Master's Thesis 'Creating Phylogenetic Networks from Clusters' available at http://resolver.tudelft.nl/uuid:4dd0bbf1-a734-4e72-af8f-de3c3651c760.

## Dependencies
The Python code dependends on the following packages:
- Biopython
- Matplotlib (only if figures are shown, i.e. if the result is saved as figure)

## How to run MSTCass
1. Download the source code from the folder MSTCass;
2. Go to the folder MSTCass;
3. From a terminal, run e.g. "python3 main.py -infile=ElusivenessFig9.nwk -outfile=result.nwk -fig_name=result".

Notes:
- The given input filename is searched in the directories '../Cass_data/restricted_grass_trees/' and '../Cass_data/restricted_grass_clusters/' automatically.
- Run "python3 main.py -h" for help
- Giving the base filename for a figure is enough; '.png' will be appended automatically.
- As path to save the figure, one may include (sub)direcories, e.g. "-fig_name=output/result"

## How to run Cass
1. Download the source code from the folder Cass_src;
2. Build the three Java files (using e.g. the command "javac *.java -d outputdir"). Note that the files were built using Java 11;
3. Go to the directory (e.g. outputdir) containing the built class files;
4. Run the algorithm, using e.g. the command "java CassAlgorithm '../Cass_data/restricted_grass_clusters/ElusivenessFig9.clu' --printretnum --timeout=5".

Notes:
- The option --printretnum prints the progress of the reticulation number found so far. The given number always indicates a lower bound. The format is 'k=K,r=R', where K is the minimum level of the output network and R is the minimum reticulation number. If a network is found, the data that are printed lastly, represent the parameters for the output network (e.g. 'k=4,r=4' for ElusivenessFig9.clu).
- The option --timeout=X ensures that the program does not run longer than X minutes.
