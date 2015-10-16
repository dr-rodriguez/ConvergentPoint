# Convergent Point Tools

Here is a variety of tools I've created to run a convergent point analysis on a group of stars.
A javascript version of this tool can be run [here](http://www.das.uchile.cl/~drodrigu/CPCalc.html).

This code is based on the methods described in de Bruijne (1999 MNRAS, 306, 381), Mamajek (2005 ApJ, 634, 1385), and Jones (1971 MNRAS, 152, 231).
The convergent points used here are listed in Rodriguez et al. (2013, ApJ 774, 101). In addition to running the 
analysis for an individual group, this tool can present the results for all stored groups simultaneously. We kindly ask that if 
this tool and its results are useful in your work you reference our paper.

The idea behind the convergent point analysis is that one takes the parameters of the object and computes the proper motions 
in directions parallel and perpendicular to the location of the convergent point. Probabilities are estimated from the magnitude 
of the perpendicular component of proper motion and associated errors (see Equation 23 in de Buijne 1999). We note that varying 
the internal dispersion and distance to the group only marginally affects the probabilities. The sigma-theta term is a measure 
of the uncertainty in the convergent point location (see Equation 6).

Results of these tool should be used with caution (see our paper for details and a comparison with BANYAN). A high probability 
does not necessarily imply the star is young and a member of the particular moving group. In particular, compare the predicted 
distance with the group distance. A large difference may suggest the target is not a real member. 
We recommended interested users compare results for their objects with the [Bayesian Analysis Tool](http://www.astro.umontreal.ca/~malo/banyan.php) developed by Malo et al. (2013).

Important Note: When comparing with BANYAN, note that this convergent point tool does not take into account measured radial 
velocities or distances whereas BANYAN can. Hence, if you wish to directly compare membership likelihoods from both tools make 
sure to use only object coordinates and proper motions.