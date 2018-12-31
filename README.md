# PromoterStrengthPredictorML
Machine Learning using Multivariate Linear Regression to Predict the Strength of a Promoter


[![DOI](https://zenodo.org/badge/105354787.svg)](https://zenodo.org/badge/latestdoi/105354787)


Quick summary

This website presents a Python based machine learning platform to predict the strength of sigma70 core promoters in Escherichia coli to ease the need for laborious and expensive experiments. Here multi-variate linear regression has been used where the parameters were optimized with gradient descent. The training data set used here is the promoter collection characterized by the Anderson lab at UC Berkeley (http://parts.igem.org/Promoters/Catalog/Anderson). The -35 and -10 motifs were extracted from RegulonDB (as described in Bharanikumar, Premkumar and Palaniappan, PromoterPredict: sequence-based modelling yields logarithmic dependence between promoter strength and sequence. PeerJ, 2018). They are available in the Datasets folder and are used by the standalone program (model.py) for predicting promoter strength.

How to Set it up:

General Dependencies

* Python Version 2.7
* Numpy
* Biopython
* Matplotlib

Dependencies to Interface with the Client

* Installation of Webpy 
* Nginx Server

Contributions:

1. model.py is the standalone. Install the latest version of numpy, Biopython and matplotplib; ensure that the path to the libraries of -35 and -10 motifs (available in the Datasets folder) in the program is correct, and execute:

    $python model.py

2. Finalp.py has the code which is fully interfaced with the web.

How to use the software/web-server:

1. Enter the -35 and -10 hexamers of the promoter in question (whose strength needs to be predicted):
2. You could data on characterized promoters to add to the dataset used for model building. First specify the number of instances you wish to add; next specify the -35 motif (string in nucleotide alphabet), -10 motif (string in nucleotide alphabet) and the promoter strength (float) for each promoter instance in that order.
3. The predicted strength is returned, along with the re-computed goodness of fit and regression surface of the possibly updated model.

Refer our manuscript for further details.

Bharanikumar, Premkumar, and Palaniappan. PromoterPredict: Sequence-based modelling of Escherichia coli sigma70 promoter strength yields logarithmic dependence between promoter strength and sequence. PeerJ, 2018

Authors:
 
* Ashok Palaniappan
* Ramit B
* Keshav Aditya R.P
