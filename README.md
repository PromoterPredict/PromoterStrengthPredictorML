# PromoterStrengthPredictorML
Machine Learning using Multivariate Linear Regression to Predict the Strength of a Promoter


[![DOI](https://zenodo.org/badge/105354787.svg)](https://zenodo.org/badge/latestdoi/105354787)


Quick summary

This website presents a Python based machine learning platform to predict the strength of sigma70 core promoters in Escherichia coli to ease the need for laborious and expensive experiments.
Here multi-variate linear regression has been used where the parameters were optimized with gradient descent. The training data set used here is the promoter collection characterized by the Anderson lab at UC Berkeley (http://parts.igem.org/Promoters/Catalog/Anderson).

How to Set it up

General Dependencies

* Python Version 2.7
* Numpy
* Biopython
* Matplotlib

Dependencies to Interface with the Client

* Installation of Webpy a Python Framework
* Nginx Server

Contribution Guidelines

* Model.py is the standalone. Install dependencies and execute at the OS prompt:

  $python model.py

* Finalp.py has the code which is fully interfaced with the web.

* HOW THE MODEL WORKS FOR PREDICTING AND ALLOWING THE USER TO ADD MORE DATASETS

* OBTAIN ALL THE INPUTS FROM THE CLIENT AS A SINGLE STRING : Eg TGTGTGTGATGA*TGATGATGTGTG#0.5&
* [THIS IS JUST A RANDOM SEQENCE TAKEN FOR THE SAKE OF THIS EXAMPLE]

* EXAMPLE OF STRING SENT FROM SERVER : TGTGTGTGATGA*TGATGATGTGTG#0.5&

* -35 AND -10 SEQUENCE THE USER ENTERS FOR WHICH HE REQUIRES THE PREDICTED STRENGTH (TGTGTGTGATGA*)
* -35 SEQUENCE : TGTGTG
* -10 SEQUENCE : TGATGA

* SEPARATOR *

* -35 AND -10 SEQUENCE THE USER ENTERS AS DYNAMIC INPUT WHICH IS TO BE ADDED TO THE DATA SET (TGATGATGTGTG#)   
* -35 SEQUENCE : TGATGA
* -10 SEQUENCE : TGTGTG
* SUPPOSE THE USER ENTERS MORE THAN 1 DYNAMIC INPUT, LETS SAY 2 INPUT SETS (TGATGATGTGTGATGATGGATTGA*)
* SEQUENCE 1 : TGATGATGTGTG
* SEQUENCE 2 : ATGATGGATTGA
* -35 SEQUENCE1 : TGATGA
* -10 SEQUENCE1 : TGTGTG
* -35 SEQUENCE2 : ATGATG
* -10 SEQUENCE2 : GATTGA

* SEPARATOR #

* THE STRENGTH FOR THE RESPECTIVE -35 AND -10 SEQUENCE THE USER ENTERS AS DYNAMIC INPUT WHICH IS TO BE ADDED TO THE DATA SET(0.5&)
* STRENGTH : 0.5
* SUPPOSE THE USER ENTERS MORE THAN 1 DYNAMIC INPUT, LETS SAY 2 INPUT SETS (0.5&0.6&)
* STRENGTH1 : 0.5
* STRENGTH2 : 0.6

Authors:
 
* Ashok Palaniappan
* Ramit B
* Keshav Aditya R.P
