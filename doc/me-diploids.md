---
title: $m_e$ in diploids
author: Arthur Zwaenepoel
geometry: margin=1in
fontsize: 12pt
---

Here is a solid and fairly transparent derivation of the formulae in
@zwaenepoel2024, without the haplodiplontic complications.

Assume a source population $A$ with known allele frequencies $q_i^\ast$
$i=1\dots L$ for variants that are locally deleterious in a sink population
$B$. We write $q_i$ for the corresponding allele frequencies in $B$.
Migration occurs from $A$ to $B$ at rate $m$.
The mean fitness among resident individuals is
\begin{align}
  \overline{W_B} 
    &= \prod_{i=1}^L (p_i^2 + 2p_iq_i(1-s_ih_i) + q_i^2(1-s_i)) \\
    &= \prod_{i=1}^L (1 - 2s_ih_i p_iq_i -s_i q_i^2) \\
    &\approx \exp\left( - \sum_i^L 2s_ih_i p_iq_i +s_i q_i^2\right)
\end{align}
The relative fitness of a migrant individual in the resident background is:
\begin{align}
  W_M &\approx \frac{1}{\overline{W_B}}\exp\left( 
          - \sum_i^L 2s_ih_i p_iq_i^\ast +s_i q_i^{\ast 2}\right)
\end{align}
write (to avoid double exponents, LaTeX...) $q_0$ for $q^\ast$, we get
factors for each locus of the form
\begin{align}
  sq^2 + 2shpq - sq_0^2 -2shp_0 q_0 
      &= sq^2 + 2shq - 2shq^2 - sq_0^2 - 2shq_0 + 2shq_0^2 \\
      &= s(q^2(1 - 2h) - q_0^2(1 - 2h) + 2h(q - q_0)) \\ 
      &= s\left[(1-2h)(q^2 - q_0^2) + 2h(q-q_0)\right] \\
      &= s\left[(1-2h)(q - q_0)(q + q_0) + 2h(q-q_0)\right] \\
      &= s(q-q_0)\left[(1-2h)(q + q_0) + 2h\right] \\
      &=-s(q_0-q)\left[(1-2h)(q + q_0) + 2h\right]
\end{align}
so we get 
  $$W_M \approx \exp\left(-\sum_i^L s(q_i^\ast-q_i)\left[(1-2h_i)(q_i + q_i^\ast) + 2h\right]\right)$$
This is also the expression in Himani's notes on heterosis.

Now consider the first generation (F1). The mean fitness among F1's is
  $$\tilde{W}_{0} \approx \exp\left(- \sum_i^L s_ih_i\left(p_iq_i^\ast + p_i^\ast q_i\right) + s_i q_i^\ast q_i\right)$$
relative to the resident, we get factors of the form
  $$sq^2 + 2shpq - sh(pq_0 + p_0q) + sq_0q = -s(q_0 - q)(h + q(1-2h))$$
So 
\begin{align}
  W_0 &\approx \exp\left(- \sum_i^L s_i(q^\ast_i - q_i)(h_i + q_i(1-2h_i))\right) \\
      &= \exp\left(- \sum_i^L s_i(q^\ast_i - q_i)(h_i(p_i-q_i) + q_i)\right)
\end{align}
Under the assumption that all migrants and their descendants cross with
residents, and that the effects of selection on allele frequencies within
these F1s, BC1s, *etc.* is negligible, than we get in the $k$th BC
generation at a generic locus
\begin{align}
  q^{(k)} = q + \frac{1}{2^{k+1}}(q^\ast - q)
\end{align}
The proportion of heterozygotes is hence
\begin{align}
  pq^{(k)} + p^{(k)}q 
&= pq^{(k)} + q - q^{(k)}q  \\
&= q^{(k)}(p - q) + q \\
&= \left(q + \frac{1}{2^{k+1}}(q^\ast - q)\right)(p - q) + q \\
&= 2pq + \frac{1}{2^{k+1}}(q^\ast - q)(p - q)
\end{align}
The proportion of homozygotes is
\begin{align}
  q^{(k)}q &= q^2 + \frac{1}{2^{k+1}}(q^\ast - q)q
\end{align}

In general for the $k$th generation backcross, we obtain
$$\tilde{W}_{k} \approx \exp\left(- \sum_i^L 
  s_ih_i\underbrace{\left(2p_iq_i + \frac{1}{2^k}(p_i - q_i)(q_i^\ast - q_i)\right)}_{\text{heterozygotes}}
  + s_i \underbrace{\left(q_i^2 + \frac{1}{2^k}q_i(q_i^\ast -
    q_i)\right)}_{\text{homozygotes}}\right)$$
so, dividing by the mean resident fitness, we get factors of the form
\begin{align}
  2shpq &+ sq^2 - 2shpq - sh(p-q)(q_0-q)/2^k -sq^2 -sq(q_0-q)/2^k \\
  &= -sh(p-q)(q_0-q)/2^k -sq(q_0-q)/2^k\\
  &= -\frac{s(q_0-q)(h(p-q) + q)}{2^k}
\end{align}
for the heteryzotes' part, and of the form
\begin{align}
  sq^2 - sq^2 - q(q_0-q)/2^k = q(q_0-q)/2^k
\end{align}
for the part coming from homozygotes. Putting everything together, this
yields:
\begin{align}
W_k
&= \exp\left(-\frac{1}{2^{k}} \sum_i^L s_i(q^\ast_i - q_i)(h_i(p_i-q_i) + q_i)\right)
\end{align}

The gff is then
\begin{align}
g &= W_M\prod_{k}^\infty W_k \\
  &= \exp\left(-\sum_i^L s_i(q_i^\ast-q_i)\left[(1-2h_i)(q_i + q_i^\ast) + 2h\right]\right)
   \exp\left(-2\sum_i^L s_i(q^\ast_i - q_i)(h_i(p_i-q_i) + q_i)\right)
   \label{eq:gff}\\
  &= \exp\left(-\sum_i^L s_i(q_i^\ast-q_i)\left[(1-2h_i)(q_i + q_i^\ast) + 2h\right]
  +2s_i(q^\ast_i - q_i)(h_i(p_i-q_i) + q_i)\right)\\
  &= \exp\left(-\sum_i^L s_i(q_i^\ast-q_i)\left[(1-2h_i)(q_i + q_i^\ast) + 2h
  +2(h_i(p_i-q_i) + q_i) \right]\right) \\
  &= \exp\left(-\sum_i^L 2s_i(q_i^\ast-q_i)\left[2h + \frac{1}{2}(3q_i - q_i^\ast)(1-2h)\right]\right)
\end{align}
in the absence of dominance, this becomes
\begin{align}
g  &= \exp\left(-\sum_i^L s(q_i^\ast-q_i)\right)
   \exp\left(-2\sum_i^L \frac{s_i}{2}(q^\ast_i - q_i)\right)\\
   &= \exp\left(-2\sum_i^L s(q_i^\ast-q_i)\right)
\end{align}
which is the same as the result for a haploid model.

When migration happens after selection, the factor $W_M$ should be dropped.
In general it appears more transparent to write the gff as in
@eq:gff.

When allele frequencies fluctuate, one can take $\Ex[g]$ as an approximation,
dropping terms that are $O(s^2)$, assuming LE, and assuming the source and
island allele frequencies are independent. The gff will then be a function
of the first two moments of the allele frequency distribution (only the
first moment in the absence of dominance).
\begin{align}
\Ex[g] &\approx \Ex[W_M]\prod_{k}^\infty \Ex[W_k] \\
  &\approx \exp\left(-\sum_i^L \Ex\left[s_i(q_i^\ast-q_i)\left[(1-2h_i)(q_i + q_i^\ast) + 2h\right]\right]\right)\\
  &\qquad\qquad \times \exp\left(-2\sum_i^L \Ex\left[s_i(q^\ast_i - q_i)(h_i(p_i-q_i) + q_i)\right]\right)\\
  &\approx \exp\left(-\sum_i^L 
    s_i\left((q^\ast_i - \Ex[q_i]) - (1-2h_i)\left[pq^\ast_i - \Ex[pq_i]\right]\right)\right) \\
  &\qquad\qquad \times \exp\left(-2\sum_i^L s_i\left(h_i(q^\ast_i - \Ex[q_i]) -
    (1-2h_i)\left[p^\ast_i\Ex[q_i] - \Ex[pq_i]\right]\right)\right)
\end{align}
This is what's in equation 10 of our genetics paper.

**Question**: how do random allele frequencies affect the assumptions
regarding the allele frequencies in the BC generations? Are we safe in
using the $\Ex[g]$ where $g$ is obtained conditional on known allele
frequencies?

Two-locus theory suggests that the effective migration rate at a neutral locus
linked to a selected one in a diploid is given by
\begin{align}
m_e 
&= m\left(1 - \frac{s(q - q^\ast)(h + (1-2h)q)}{r + m + s(h + (1-2h)q)(p-q)}\right) \\
&= m\left(1 - \frac{s(q - q^\ast)(h(p-q) + q)}{r + m + s(h(p-q) + q)(p-q)}\right)\\
&= m\exp\left(- \frac{s(q - q^\ast)(h(p-q) + q)}{r + m + s(h(p-q) + q)(p-q)}\right)
\end{align}
For the case where $r \approx 0.5 \gg m,s$ this yields the same prediction
as above (barring the contribution from the diploid migrant).

If we take the exponent in the last line, and expand in powers of $s$, we
find 
  $$-\frac{s}{m + r} \left(2 h q^{2} - 2 h q q_{0} - h q + h q_{0} - q^{2} + q q_{0}\right) + O\left(s^{2}\right)$$
Taking expectations and rearranging, we obtain to first order in $s$
  $$-\frac{s}{m + r} \left(h\left(q^\ast - \Ex[q]\right) - (1-2h)\left[p^\ast \Ex[q] - \Ex[pq]\right]\right)$$
This suggests a gff approximation of the form
\begin{align}
\Ex[g(x)] &\approx \Ex[W_M]\prod_{k}^\infty \Ex[W_k] \\
  &\approx \exp\left(-\sum_i^L 
    s_i\left((q^\ast_i - \Ex[q_i]) - (1-2h_i)\left[pq^\ast_i - \Ex[pq_i]\right]\right)\right) \\
  &\qquad\qquad \times \exp\left(-\sum_i^L\frac{s_i\left(h_i(q^\ast_i - \Ex[q_i]) -
    (1-2h_i)\left[p^\ast_i\Ex[q_i] - \Ex[pq_i]\right]\right)}{m + r(x,x_i)}\right)
\end{align}

However, this ignores the fact that we often will have $m=O(s)$ and
$r=O(s)$ as $s\rightarrow 0$. Assuming $m=O(s)$ and $r=O(s)$, we get
  $$-s\left(\frac{q q_{0} - q^{2} + h q_{0} - h q - 2 h q q_{0} + 2 h q^{2}}{
m + r + s \left(4 h q^{2} - 4 h q + h - 2 q^{2} + q\right)}\right) + O\left(s^{2}\right)$$
Taking expectations is now more tricky because of the $q$'s in the
denominator. Proceeding naively with the numerator and denominator
separately, we get
  $$-\frac{s\left(h\left(q^\ast - \Ex[q]\right) - (1-2h)\left[p^\ast \Ex[q] - \Ex[pq]\right]\right)}{m + r + s(h - \Ex[q] + (1-2h)2\Ex[pq])}$$
This suggests a gff approximation of the form
\begin{align}
\Ex[g(x)] &\approx \Ex[W_M]\prod_{k}^\infty \Ex[W_k] \\
  &\approx \exp\left(-\sum_i^L 
    s_i\left((q^\ast_i - \Ex[q_i]) - (1-2h_i)\left[pq^\ast_i - \Ex[pq_i]\right]\right)\right) \\
  &\qquad\qquad \times \exp\left(-\sum_i^L\frac{s_i\left(h_i(q^\ast_i - \Ex[q_i]) -
    (1-2h_i)\left[p^\ast_i\Ex[q_i] - \Ex[pq_i]\right]\right)}{
    m + r(x,x_i) + s_i\left(h_i - \Ex[q_i] + (1-2h_i)2\Ex[pq_i]\right)}\right)
\end{align}
which is equation 5 in the Genetics paper.

Consider $h=1/2$ to obtain a haploid/no dominance model:
\begin{align}
E[g(x)] &\approx E[W_M]\prod_{k}^\infty E[W_k] \\
  &\approx \exp\left(-\sum_i^L s_i(q^\ast_i - E[q_i])\right) 
   \exp\left(-\sum_i^L\frac{s_ih_i(q^\ast_i - E[q_i])}{
    m + r(x,x_i) + s_ih_i\left(1 - 2E[q_i]\right)}\right)\\
  &\qquad \text{where } h_i = 1/2
\end{align}

What does this entail for a single selected haploid locus at equilibrium
frequency $m/s$, fixed in the mainland?
\begin{align}
E[g(x)] 
  &\approx \exp\left(-s(1 - m/s)\right) \exp\left(-\frac{s(1 - m/s)}{
    m + r + s\left(1 - 2m/s\right)}\right)\\
  &\approx \exp\left(-(s - m)\right) \exp\left(-\frac{s - m}{r + s - m}\right)\\
  &\approx \exp\left(-(s - m)\left(1 +\frac{1}{r + s - m}\right)\right)\\
  &\approx \exp\left(\frac{(m - s)(1-m+r+s)}{r + s - m}\right)\\
  &\qquad\text{where } m < s
\end{align}
as $m \rightarrow 0$, this becomes
$$\exp\left(\frac{-s(1-r-s)}{r+s}\right) \approx
\exp\left(\frac{-s}{r+s}\right) \approx 1 - \frac{s}{r+s} = \frac{r}{r+s}$$
Which is the result of @petry1983.


## References
