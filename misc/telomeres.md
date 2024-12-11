## Considering telomeres and centromeres

#### IBD segment detection
---

`hap-ibd.jar` and `ibd-ends.jar` are not good at detecting IBD segments near the ends of the genetic map. The pipeline tries to address this issue by via the parameter `TELOSIGMA`.
- `TELOSIGMA` is an initial removal of outlier regions before computing a mean and standard deviation. It removes regions with such and such standard deviations above the median.
    - The IBD segment detection methods are also not good near centromeres. The initial removal with `TELOSIGMA` may address this issue
    - The initial removal with `TELOSIGMA` may also address excess variance due to regions that are strongly selected for. In European humans, the LCT signal is enormous (likely genome-wide significant), and without the initial removal can create excess variance without the initial removal.
- The genetic maps may be poor estimates as well near telomeres and centromeres.
- I also try to address this by triming the ends of genetic maps by the size of the parameter 'SCANCUTOFF'.
    - There would otherwise be an issue with truncation when detecting segment lengths.

#### Small chromosomes
---

I recommend against analyzing chromosomes of length less than 10 cM.
- The pipeline assumes that chromosomes are consecutive ordered from `CHRLOW` to `CHRHIGH`.
- If after ignoring chromosomes smaller than 10 cM, resulting in a gap in the chromosome ordering, do the following.
    - Create file pointers `ln -s old_file new_file` for genetic maps and VCFs
    - Modify the *.yaml file accordingly.