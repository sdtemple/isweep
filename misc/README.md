The *.jar files are in-house apps from https://github.com/browning-lab.
- These are not being supported w/ documentation or feedback.
- They do simple manipulations of text data.
    - filter-lines.jar is important to the analysis pipeline.
        - It subsets the text file based on values in a single column.
        - Consider it a newer version of https://faculty.washington.edu/browning/beagle_utilities/utilities.html
    - add-uniform-err.jar and remove-phase.jar are important to simulation study.
        - add-uniform-err.jar changes some alleles (in vcf) at a user-defined rate.
        - remove-phase.jar changes say 0|1 or 1|0 to 0/1