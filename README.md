# drugCIPHER
This repository contains the implementation of drugCIPHER, 
a linear regression framework to predict drug-target relations.
For more details please see 
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011764

## Introduction
_drugCIPHER_ is a linear regression framework to integrate heterogenous drug 
similarities with protein interation network data to acurately predict drug-target 
relations. Three linear regression models are proposed respectively using drug 
therapeutic similarity, chemical similarity and their combination as responses and 
network distance as predictors.


## Code
The _code_ folder contains four matlab files for drugCIPHER. 
* _drugCIPHER_SingleS_Validation.m_ performs leave-one-out cross-validation on drug target prediction with single drug similarity matrix as input
* _drugCIPHER_SingleS_Overall.m_ performs drug target overall prediction with single similarity matrix 
 for all input drugs. Candidate targets are treated as all genes/proteins in the protein interation network 
* _drugCIPHER_MS_Validation.m_ performs leave-one-out cross-validation with two similarity matrices
*  _drugCIPHER_MS_Overall.m_ performs drug target overall prediction with two drug similarity matrices
