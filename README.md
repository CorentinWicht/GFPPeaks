# GFPPeaks

The MATLAB script in this repository enable the semi-automatic identification of [Global Field Power (GFP)](https://www.hindawi.com/journals/cin/2011/813870/) peaks on Event-Related Potentials (ERP) files. 

The method is based on: [Wicht, C.A., De Pretto, M. & Spierer L.; Neural correlates of expectations-induced effects of caffeine intake on executive functions; Cortex Registered Reports, *stage 1 in principle acceptance*](https://osf.io/sudnm/)


**⚠️ OF NOTE: The analysis script can currently only import .ep EEG files (see [Cartool](https://sites.google.com/site/cartoolcommunity/)).**


## Cite the repository
C.A. Wicht, GFPPeaks, (2020), GitHub repository, https://github.com/CorentinWicht/GFPPeaks \
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4032759.svg)](https://doi.org/10.5281/zenodo.4032759)

## Table of Contents
- [GFPPeaks](#gfppeaks)
  * [Cite the repository](#cite-the-repository)
  * [Getting Started](#getting-started)
    + [1. Files Extensions and Processing Mode](#1-files-extensions-and-processing-mode)
    + [2. Loading Files](#2-loading-files)
    + [3. Design Definition](#3-design-definition)
    + [4. Unblinding](#3-unblinding)
    + [5. Data Files Selection](#4-data-files-selection)
    + [6. Include Files](#5-include-files)
    + [7. ERP Parameters](#6-erp-parameters)
    + [8. ERP Components Definition](#7-erp-components-definition)
    + [9. GFPPeaks Windows](#8-gfppeaks-windows)
      - [9.1 Main Window](#81-main-window)
      - [9.2 Data Table](#82-data-table)
    + [10. Exported Figure](#9-exported-figure)
  * [Author](#author)
  * [License](#license)
  * [Acknowledgements](#acknowledgements)
  * [Fundings](#fundings)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



## Getting Started

When on the [GFPPeaks startup page](https://github.com/CorentinWicht/GFPPeaks), start by clicking on `⬇️ Code` on the top right of the screen and then `Download ZIP` to download the whole repository (alternatively you can also clone it). 
Then, you will need to run the following script: ```GFPPeaks.m```


You will find below a step by step explanation on how to run each script in MATLAB (by clicking on the `▶️ Run` button or on your keyboard's button `⌨️ F5`).


### 1. Files Extensions and Processing Mode

![](tools/screenshots/Settings.png)

You can define first define the extension of the files you would like to load (currently .ep or .eph). 

Then, you can choose between a : 
1) semi-automated 
2) manual processing modes 
3) or loading previous data \
The semi-automated method looks for the GFP peak inside each component's upper and lower bounds as defined in [Chapter 7](#7-erp-components-definition), *see* [Chapter 8.1](#81-main-window) *for an example of the main GFP Peaks window.* \
The manual method requires you to enter manually the GFP peaks for each files (identified for e.g. with CARTOOL) in the table detailed in [Chapter 8.2](#82-data-table). 
Loading data will require you to have ran the script previously to get the following .mat file(s):

![](tools/screenshots/LoadMAT.png)


### 2. Loading Files

![](tools/screenshots/LoadFiles.png)

Here you are asked to provide the most upper folder containing all your ERP files with .ep extension.


### 3. Design Definition

![](tools/screenshots/Design.png)

The third prompt enables the definition of the statistical design.

These are the designs that are currently accepted:
```
* 1 between-subjects factor
* 1 within-subjects factor
* 1 between- & 1 within-subjects factors
* 2 between-subjects factors
```
There can only be **up to 3 levels for each factor !**

There is already an example that shows you how to fill the first two columns, namely on the... :
```
* First line you should indicate whether Factor1 is a Between-subjects (B) or Within-subjects (W) factor
* Second line you should give a name to the Factor
* Third and + lines you need to define the name of each levels.\
```
**This is case-sensitive since the name you pick as levels should match either a pattern in your EEG files or the name of subfolders in which you stored the files.**

Additionally, you have the possibility to ignore specific `folder(s)` and/or `file(s)` that are inside the EEG files-containing upper folder.\
For example, if you would like to ignore every EEG files that contains the pattern `_EXCLUDE`, just write this pattern in the first box under the `IgnoreFiles` header.

**Once you are done, press on the cross `X` on the top-right corner and the code will resume**

### 4. Unblinding

![](tools/screenshots/Unblinding.png)

Since the individual GFP peaks will be averaged for each factor's level, you can here unblind the conditions assignment.


![](tools/screenshots/AllIndivUnblind.png)

You further need to determine whether you  want to replace the factor's levels names by:
1) their unblinded names for all files (All):

![](tools/screenshots/AllUnblind.png)

2) or each individual file need to be renamed separately (Individual, ! requires Excel !):

![](tools/screenshots/ExcelUnblind.png)

**CAREFUL NOT TO CLOSE THE EXCEL FILE, JUST PRESS "OK" WHEN DONE**


### 5. Data Files Selection

![](tools/screenshots/RestrictData.png)

Here you have the possibility whether you would like to process all the data in the folder path you provided in [Chapter 1](#1-loading-files) or only a subset of the files. \


### 6. Include Files

![](tools/screenshots/IncludeFiles.png)

Here you have the possibility to decide which data you like to process. \
Close the window once you are satisfied with your choice.

![](tools/screenshots/SucessfullyLoaded.png)

In MATLAB's command window you will see the list of correctly loaded files according to your selection in the previous prompt.


### 7. ERP Parameters

![](tools/screenshots/ERPParams.png)

In this prompt you need to:
```
1) define the epoching interval of your files (e.g. -100 to 700ms), in MILLISECONDS (!). 
2) report the sampling rate (e.g. 1024Hz).
3) Define a folder name where to save the results (e.g. "Output").
```

### 8. ERP Components Definition

![](tools/screenshots/ComponentsDef.png)

This is one of the most important prompt ! \
On each of the four lines you have the possibility to define a specific ERP component of interest. \
Each column should be filled accordingly:
```
Column 1: Type in the name of the component (e.g. N2)
Column 2: Define the lower bound in ms for the component of interest (e.g. 200ms)
Column 3: define the upper bound in ms for the component of interest (e.g. 350ms)
```

These informations will be used in the GFP peaks main display in [Chapter 8.1](#81-main-window)


### 9. GFPPeaks Windows
#### 9.1 Main Window

![](tools/screenshots/MainWindow.png)

This is the main window which will enable you to determine whether the preselected GFP peak in the data table (*see* [Chapter 8.2](#82-data-table)) is correct. \
The window displays *from left to right, top to bottom*:
```
1) The ERP including all channels (each channel is represented as one line) [top left].
2) The GFP activity for epoch of interest [bottom left].
3) The topographical plot (topoplot) [right].
```
The vertical line represent, for all graphs, the position of the maximum GFP identified by the script. \
Next to the vertical line, you have the exact peak position in ms and in time-frames (TF; i.e. dependent of your sampling rate). \
The vertical blue bar represents the component of interest's time range. 

On the far right you will find the topography at the GFP peak (indicated in TF in the title). \
You have the possibility to use the slide bar below to inspect the topography of the whole ERP period while this will also adjust the vertical line in all other graphs.


#### 9.2 Data Table

![](tools/screenshots/GFPPeakData.png)

This table will popup next the main figure from [Chapter 8.1](#81-main-window). \
**⚠️ DO NOT CLOSE THIS TABLE BEFORE YOU ARE FINISHED WITH THE CURRENT FILE ⚠️**.

The GFP peak for the current file will be prefilled if you selected the semi-automatic method while the table will be completely empty if you selected the manual method (*see* [Chapter 2](#2-processing-mode)).

For the semi-automatic method, either you leave the prefilled value if you think it is correct or you can change it manually. \

**Once you are satisfied with your results you can close the data table and the results will be saved**.

### 10. Exported Figure

![](tools/screenshots/Results.png)

You are now done with the script and based on what you gave in [Chapter 6](#6-erp-parameters) as folder name for the outputs, you fill find the figure above.\
This figure shows the mean GFP activity for the provided ERP epoch, separated by conditions/groups (i.e. depending on your design). \
The vertical red line represents the peak (i.e. maximum) of the mean of each individual GFP, while the red rectangle represents +- 1 Standard Deviation (SD) over the each individual GFP.\
Finally, below the GFP graph you will find the topography averaged over the period included in +- 1 Standard Deviation (SD) around the mean peak (i.e. time range of averaging indicated in TF in the title).

## Author
[**Corentin Wicht**](https://www.researchgate.net/profile/Wicht_Corentin)\
*SNSF Doc.CH PhD student*\
*corentin.wicht@unifr.ch, corentinw.lcns@gmail.com*\
*[Laboratory for Neurorehabilitation Science](https://www3.unifr.ch/med/spierer/en/)*\
*University of Fribourg, Switzerland*

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

See the [LICENSE.md](LICENSE.md) file for details

## Acknowledgements
PD Dr. Lucas Spierer, Director of the [Laboratory for Neurorehabilitation Science (LNS), Section of Medicine, Faculty of Science and Medicine, University of Fribourg, Switzerland](https://www3.unifr.ch/med/spierer/en/) provided substantial support and advices regarding theoretical conceptualization as well as access to the workplace and the infrastructure required to successfully complete the project.

## Fundings
This project was supported by [Swiss National Science Foundation](http://www.snf.ch/fr/Pages/default.aspx) grants:
* [#P0LAP1_181689](http://p3.snf.ch/project-181689) to Corentin Wicht
* [#320030_175469](http://p3.snf.ch/project-175469) to PD Dr. Lucas Spierer
