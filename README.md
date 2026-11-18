*This repo contains data and R code that is sufficient to reproduce the analysis in:*

**Aging Patterns in the Blood Metabolome of Roe Deer in the Wild**<br>
by *Benjamin R Harrison, Benjamin Rey, Lu Wang, Richard Johnson, Michael J MacCoss, François Débias, Revelli Paul, Maryline Pellerin, Daniel. E. L. Promislow, Jean-Michel Gaillard, and Jean-François Lemaître*<br>
2026, *Submitted to Aging Cell* 

The R code files (**in bold** below) will use ether data from the *data* folder, or data files denerated by the other R files.  
In order to complete the analysis, you will need to run some of the R files before the others.  The data files used in each R code file are listed, alond with their source; either the data folder, or 'internal', meaning they are made by another R file.  Any data that the R files generate that are used by the other files are also listed. <br>

**data merge and normalization.R**<br>
uses: <br>
roeDeerAreaReport3-METABOLOME.csv [data folder] <br>
RoeDeer-LIFE HISTORY DATA-V6.xlsx [data folder] <br>

generates: <br>
dataforUnivariateModels
<br>

**PCA.R**<br>
uses: 
dataforUnivariateModels [internal]<br>

does not generate files for other code<br>

**univariateAnalysis.R**<br>
uses: <br>
dataforUnivariateModels [internal]<br>

does not generate files for other code<br>

**RoeDeerClock.R**<br>
uses: <br>
clockBuildingFucntion [data folder]<br>
RoeDeer-LIFE HISTORY DATA-V6.xlsx [data folder]<br>
dataforUnivariateModels [internal]<br>

generates: <br>
clockTraining
<br>

**RoeDeerClock_Sex_Specific.R**<br>
uses: <br>
clockBuildingFucntion [data folder]<br>
RoeDeer-LIFE HISTORY DATA-V6.xlsx [data folder]<br>
dataforUnivariateModels [internal]<br>

does not generate files for other code<br>

**longitudinalAnalysis.R**<br>
uses: <br>
RoeDeer-LIFE HISTORY DATA-V6.xlsx [data folder]<br>
dataforUnivariateModels [internal]<br>

does not generate files for other code<br>


