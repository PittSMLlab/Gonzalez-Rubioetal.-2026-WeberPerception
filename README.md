# Weber's Law in walking: sensory scaling is observed in multi-sensory, dynamic tasks. 
Gonzalez-Rubio, M., Iturralde, P.A. & Torres-Oviedo, G. Weber’s Law in walking: sensory scaling is observed in multi-sensory, dynamic tasks. Sci Rep (2026). https://doi.org/10.1038/s41598-026-54948-5

## Data Organization
Data from the experiment are made available in the `Data` folder of this repository. There are two types of files:
- `Data.csv` stores a table in CSV format containing task information information for all 39 participants,
including the presented belt speed difference (ΔV, in mm/s), binary choice response (left or right leg
perceived as slower), reaction time (time from auditory start cue to button press), participant ID, block
number, and speed condition.
- `JND_Choices.mat` contains the just noticeable differences (JNDs) estimated from participants' choice
data, derived from the generalized linear mixed-effects model (GLMM) described in the manuscript.

A subfolder within `Data` contains the JNDs estimated from reaction time data alone for each group
(Slow, Comfortable, and Fast), computed via the drift-diffusion model (DDM). These files are used in
conjunction with the Python code described below.

Original raw data files (.mat) will be made available upon request.

## Analysis
Figures 2 and 3 can be reproduced using **Weber_Figures.m**:
- Make sure your "CurrentFolder" in MATLAB is set to where the `Data` folder is located.
- Helper functions `Boxplot.m` and `probeColorMap.m` must be in the same folder as the main script
or added to the MATLAB path.
- MATLAB has a web-based access option, however a license and account are required.

Figure 4 can be reproduced using **DDM_RT_fit.ipynb** and **utils.py**:
- No license is needed to run the Jupyter Notebook.
- `utils.py` must be in the same directory as the notebook.
- Data for this analysis can be found in the `Data` folder as `Data.csv`.
