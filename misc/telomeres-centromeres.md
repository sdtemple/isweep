## IBD segment detection accuracy near telomeres and centromeres

### Telomeres
---

`hap-ibd.jar` and `ibd-ends.jar` are not good at detecting IBD segments near the ends of the genetic map. The pipeline tries to address this issue by via the parameter `TELOSIGMA`.
- `TELOSIGMA` is an initial removal of outlier regions before computing a mean and standard deviation. It removes regions with such and such standard deviations above the median.
    - The IBD segment detection methods are also not good near centromeres. The initial removal with `TELOSIGMA` may address this issue
    - The initial removal with `TELOSIGMA` may also address excess variance due to regions that are strongly selected for. In European humans, the LCT signal is enormous (likely genome-wide significant), and without the initial removal can create excess variance without the initial removal.
- The genetic maps may be poor estimates as well near telomeres and centromeres.
- I also try to address this by triming the ends of genetic maps by the size of the parameter 'SCANCUTOFF'.
    - There would otherwise be an issue with truncation when detecting segment lengths.

### Centromeres
---

In some human autosomes, there are large centromeric gaps in the GRCh38 reference genome. These gaps in cM and bp can result in inaccurate IBD segment detection. The 'TELOSIGMA' initial removal step is supposed to partially address this. The concern about IBD rates near centromeres may be less of an issue using a map without large gaps. We have not investigated the accuracy of IBD segment detection in highly repetitive regions.