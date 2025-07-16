
Two neat "genetic architecture hacks":

## Effective migration rates with complete divergence

For $L$ divergently selected loci with selection and dominance
coefficients $s_i$ and $h_i$, $i=1, \dots,L$, where the two
populations are near complete fixation for their respective locally
beneficial allele, we get, at position $x$ in the genome, the
following expression for $m_e$:
  $$m_e(x) = m \times \exp\left( -\sum_{i=1}^L \frac{s_i h_i}{m + r(x,x_i) + s_ih_i} \right)$$
This will be a reasonable approximation as long as $N_e sh \gg m$ and
$N_e sh \gg 1$, so that selection is much stronger than both drift and
migration.

This can be obtained by considering the reproductive value of a
migrant individual in the resident background, or by an *ad hoc*
extrapolation of two-locus theory to the multi-locus setting, see
@barton1986, @aeschbacher2017, @sachdeva2022, @zwaenepoel2024.

## Effective population size under BGS

In the "background selection regime" (see @good2014genetic), we can
use the following expression [@hudson1994gene]:
  $$N_e(x) = N \times \exp \left( -\sum_{i=1}^G \frac{u_i}{
      s_i\left(1 + r(x,x_i)\frac{(1 - s_i)}{s_i}\right)^2} \right)$$
Where there are $G$ sites with recurrent deleterious mutations at rate
$u_i$ and (haploid) selection coefficient $s_i$, $i=1, \dots, G$.
This should hold more or less whenever it provides reasonable
$N_e$'s -- see @good2014genetic, @buffalo2024.

This can be derived, more or less, from the structured coalescent, where the
different subpopulations are deleterious backgrounds, see
@hudson1994gene, @charlesworth1997, @nordborg1997structured, @charlesworth2010,
@good2014genetic, @santiago1995, @santiago1998, @santiago2016, @buffalo2024.

## How well do they work?

@zwaenepoel2024 showed that the $m_e$ theory can predict equilibrium
allele frequencies at selected loci very accurately when linkage
between selected loci is not too tight.
To what extent the theory allows to predict genome-wide *neutral*
divergence/differentiation was not studied in depth. In particular,
predictions very close to barrier loci are note xpected to be
accurate.

In addition, it is not clear how BGS and divergent selection interact.
In the regime where the Hudson & Kaplan theory applies with respect to
BGS, can we account for both processes by rescaling $N$ and $m$ in the
neutral theory? This remains to be studied in depth.

It would be interesting to study these questions systematically,
focusing on the regimes with near-complete divergence and HK-like BGS.

![
Predictions based on $m_e$ and BGS theory for forward simulations
using SLiM.
](/home/arzwa/vimwiki/build/img/2025-06-19.pdf)


## References
