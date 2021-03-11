# README #

Conditioned Latin hypercube sampling (cLHS) was proposed by [Minansy and McBratney (2006)](http://www.sciencedirect.com/science/article/pii/S009830040500292X). The aim of cLHS is to geographically locate samples (it was intended for soil sampling but can be applied to other sampling contexts too) such that the empirical distribution functions of the ancillary information associated with the samples are replicated, with a constraint that each k-tuple of ancillary information has to occur in the real world.

When using cLHS one often asks the question of how many samples should be collected? Other times, particularly when in the field and finding  out that a sample location selected by cLHS can not be visited because of any number of reasons; they may ask, where else should i go then? Ideally in this situation the person would want to sample an alternative site that is very similar to the site (based on the multivariate information) that was unable to be visited. Another situation is when planning a new survey and one wants to use cLHS, but wants to be able to take into consideration prior or existing samples that have been collected from their sampling domain. Taking the existing samples into consideration may lessen the work required in conducting the new survey or better target areas that have been under-sampled.  


### What is this repository for? ###

This repository contains a few workflows all developed in R, that answers those particular questions described above. The workflows exploit the use of the [clhs](https://cran.r-project.org/web/packages/clhs/index.html) R package developed by Pierre Roudier (Landcare Research New Zealand).


### Who do I talk to? ###

* Brendan Malone (brendan.malone@sydney.edu.au)