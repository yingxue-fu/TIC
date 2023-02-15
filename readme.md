# TMT interference correction using TMTc ions 
## Download the program
Download all files by clicking "Code" and then "Download ZIP". Unzip the file to a directory and change the working directory to it.

`cd path/to/the/directory`
## Set up conda environment
1. Follow the instructions (https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) to install conda in your operating system.
2. Use the "environment.yml" file to create the conda environment.
    
    `conda env create --file environment.yml`
3. Activate the conda environment.
    
    `conda activate tmtcNOISEcorrc`

## Run the program
There are three seperate steps for running the program, i.e., 1) TMTpro reporter ion-based quantification, 2) TMTproC ion-based quantification, and 3) TMTpro reporter ion noise correction. For each step, there is a paramter file (located in the "Programs" folder) used to set up the input and output file path and other required parameters. 

To test the program using the example input files, you need to download the corresponding mzXML file from here: https://drive.google.com/file/d/1ftBng_v8OpNbBdv7vydeVIL4xrBGuuf9/view?usp=sharing, and put it in the \Example_Input folder.

Then run the followling commands step by step.

1. TMTpro reporter ion-based quantification

    `python Programs\01_Reporter_based_quan\jump_q_reporter.py Programs\jump_q_01_reporter.params`

2. TMTproC ion-based quantification

    `python Programs\02_TMTc_based_quan\jump_q_tmtc.py Programs\jump_q_02_tmtc.params`

3. TMTpro reporter ion noise correction

    `python Programs\03_Reporter_correction\jump_q_noise_corrc.py Programs\jump_q_03_noise_corrc.params`

After the program finished running, you can check the results in the "test_Output" folder.

Note: The formats of the identification results, i.e., tables contain peptide-spectrum matches (PSMs), peptides and proteins, are based on output files from JUMP suite. You may need to adjust your input files accordingly if you want to use the identification results from other software. 


