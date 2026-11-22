# Exponential variance of paired determinant estimates for optical circuit readout

Daniel Fredriksen, Quantyra Inc. · [Manuscript PDF](manuscript.pdf) · [LaTeX source](manuscript.tex) · [ORCID](https://orcid.org/0009-0002-9123-7435)

**Cite version 1.0.0.** [Zenodo DOI 10.5281/zenodo.24430932](https://doi.org/10.5281/zenodo.24430932) · [All versions](https://doi.org/10.5281/zenodo.24430931)

This preprint gives a connected family of optical circuits on which an unbiased, shared-sign determinant estimate of a one-bit acceptance probability has exponential variance. The bound holds for every positive entry weighting fixed before sampling. Deleting any one central controlled-sign gate reverses the promised decision. A companion family uses independent control preparations. Both families admit efficient classical decision algorithms.

## Exact checks

Run from the repository root with Python 3.10 or later and assertions enabled. Both scripts use the standard library.

~~~sh
python -B checks/paired_moments.py
python -B checks/live_cycles.py
~~~

The first program contracts the full four-determinant second moment exactly for both basis inputs of the three-qubit example. Their acceptance probabilities are 25/169 and 144/169, and their normalized variances are approximately 419.484581 and 363.480837. The second builds the complete optical matrices for one triple of gates in each growing family, including preparation and reversal. It checks connectedness, all eight complementary derivatives in the arm average, and the cycle inequality at asymmetric rational weights. It also builds the 78-port circuit with three triples and verifies the products and boundaries of all 512 masks. The all-size lower bound is proved in the manuscript.

## Rebuild the PDF

Use Tectonic 0.17.0. Set the environment variable SOURCE_DATE_EPOCH to 1795219200 for the distributed PDF, then run the command below.

~~~sh
tectonic manuscript.tex
~~~

The manuscript and accompanying materials are licensed under [CC BY 4.0](LICENSE). The [review disposition](review.md) describes the automated reviews. The [Lean certification plan](lean-certification-plan.md) records prospective proof obligations.
