## Different average recombination rates

The haplotype-based methods in `workflow/model-selection` use base pair parameter settings in light of a 1.0 cM $\approx$ 1Mb conversion.

You may want to rescale the following parameters if the average recombination rate of your species is very different from that of humans.

`workflow/scan-selection`
- FIXED:ISWEEP:MBBUF

`workflow/model-selection`
- FIXED:ISWEEP:WINSIZE
- FIXED:ISWEEP:WINSTEP
- FIXED:ISWEEP:HAPSIZE
- FIXED:ISWEEP:HAPSTEP 

