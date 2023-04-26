Description, verbal summary of simulation studies

I perturb 1 parameter at a time
Everything else is fixed
Default: s=0.03, p0=0.5, Ne=bottleneck, ongoing sweep, mult inheritance, de novo sweep, 3.0 cM

SelectionCoefficient
  Studying different selection coefficients
    0.02,0.03,0.04
  Empirically unbiased
  Adequate coverage
  Interval widths
    Smallest for 0.03
    Largest for 0.04
    Okay for 0.02

TrueNe
  Simulated, inferred Ne is TrueNe; different Ne
    Bottleneck, constant 100k, exponential growth
  Performed usual inference
  Empirically unbiased
  Adequate coverage
  Interval widths
    Tight

StrongNe
  True Ne for coalescent + IBD simulation is bottleneck
  Performed inference with misspecified Ne
    Exponential growth
    Constant 100k
  Horribly biased
  Inadequate coverage

ModerateNe
  True Ne for coalescent + IBD simulations in bottleneck
  Performed inference with Ne scaled uniformly by a constant
    0.8,0.9,1.0,1.1,1.2
  Biases in predictable direction
  Inadequate coverage

StandingVariation
  Selection from standing variation
    0.1,0.05,0.025,0.01,0
  Performed inference assuming de novo sweep
  Empirically unbiased
  Adequate coverage
  Interval widths
    Increases as standing variation increases

StandingVariationV2
  Selection from standing variation
    0.1,0.05,0.025,0.01,0
  Current allele frequency p(0) is 0.25 instead of 0.5
  Performed inference assuming de novo sweep

NotOngoingSweep
  Backward time at which sweep ended (in generations)
    0,10,20,50,100
  Performed inference assuming de novo ongoing sweep
  Empirically unbiased
  Adequate coverage
    Except for tau=100 setting
    This exceptional behavior -> IBD detectable for recent 100 gens
  Interval widths
    Increases as standing variation increases

AlleleFreq
  Different p(0)
  Performed inference with true p(0)
  Empirically unbiased
  Issue in scripting for coverage
  Interval width
    Decreasing as p0 increases

WrongAlleleFreq
  Different p(0)
  Performed inference with p(0) shifted +/-
    -0.2,-0.1,-0.05,-0.02,0.00,0.02,0.05,0.10,0.20
  Empirically unbiased
  Issue in scripting for coverage

Inheritance
  Multiplicative, additive, dominance, or recessive
  Performed inference with true inheritance model
  Empirically unbiased
    Except for recessive case
    Is this a scripting error?
  Adequate coverage
    Except for recessive case
  Interval widths
    Smallest for multiplicative, additive
    Huge for recessive
    Okay for dominant

MisInheritance
  True inheritance is dominance
  Performed inference with other inheritance models
  Moderate bias for multiplicative, additive
    Huge bias for recessive
  Inadequate coverage

CentiMorgan
  Different cM threshold for IBD
    1.0,2.0,3.0,4.0
  Performed usual inference
  Empirically unbiased
    Except for 4.0; not enough IBD segments?
  Adequate coverage
    Except for 1.0; too much noise?
  Interval widths
    Smallest for 2.0, 3.0
