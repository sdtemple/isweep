API
===

Here are the main functions a user might consider. I intend the users to import all functions with ``from isweep import *``. I ordered the functions according to their importance in estimating selection coefficients. Some of the parameters with types ``numpy.array`` can accept type ``list``. There are some non-essential functions not listed; for example, functions in the ``isweep.utilities`` module to create complex `*.ne` files. "Behind the hood", object-oriented programming (Class objects) and dynamic programming help simulate the IBD of haplotypes with a fast algorithm. 

.. autofunction:: isweep.slow.empty_function

.. autofunction:: isweep.inference.read_ibd_file

.. autofunction:: isweep.utilities.bin_ibd_segments

.. autofunction:: isweep.utilities.read_Ne

.. autofunction:: isweep.inference.chi2_isweep

.. autofunction:: isweep.inference.chi2_labeled_isweep

.. autofunction:: isweep.coalescent.simulate_ibd

.. autofunction:: isweep.coalescent.simulate_ibd_constant

.. autofunction:: isweep.coalescent.simulate_ibd_isweep

.. autofunction:: isweep.coalescent.simulate_ibd_isweep_tv

.. autofunction:: isweep.slow.simulate_ibd_wf

.. autofunction:: isweep.slow.simulate_ibd_isweep_wf

.. autofunction:: isweep.coalescent.wright_fisher

.. autofunction:: isweep.coalescent.basic_coalescent

.. autofunction:: isweep.coalescent.varying_Ne_coalescent

.. autofunction:: isweep.coalescent.walk_variant_backward

.. autofunction:: isweep.coalescent.walk_variant_forward

.. autofunction:: isweep.coalescent.walk_variant_backward_tv

.. autofunction:: isweep.inference.bootstrap_standard

.. autofunction:: isweep.inference.bootstrap_standard_bc

.. autofunction:: isweep.inference.bootstrap_percentile

.. autofunction:: isweep.inference.when_freq

.. autofunction:: isweep.inference.bootstrap_freq

.. autofunction:: isweep.outgroups.make_ibd_graph

.. autofunction:: isweep.outgroups.diameter_communities

.. autofunction:: isweep.utilities.write_Ne

.. autofunction:: isweep.utilities.make_constant_Ne

.. autofunction:: isweep.utilities.make_exponential_Ne

.. autofunction:: isweep.coalescent.probability_ibd

.. autofunction:: isweep.coalescent.probability_ibd_isweep

.. autofunction:: isweep.utilities.big_format_distribution
