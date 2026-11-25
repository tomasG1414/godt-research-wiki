---
type: note
projects: [godt]
created: 2026-11-25
updated: 2026-11-25
status: current
---

# Real samples on the $a\,\delta n/\lambda$ axis

Places real samples on the axis that [[multiple-scattering-check]] indexes its results
by, to answer the open item on the [[godt]] hub: how much of the multi-shear design
study applies to actual specimens, including the two that [[s42005-025-02416-3]]
imaged.

**Caveat up front.** Apart from [[ziemczonok2019]], [[chen1998]], [[s42005-025-02416-3]]
and [[2024-06-14-598984-polarization-hilbert-holotomography]], the literature below is
**web-sourced and not ingested**. §3 says for each paper whether the abstract or the
full text was reached and by which route. No number here comes from memory. Sizes and
medium indices that a paper did not give are marked **assumed** in the table and are
*(inference)*. The $x$ column is arithmetic on those inputs (`lab/realsample_table.py`);
nothing was simulated for this note.

## Headline

1. **The axis has two readings that differ by a factor 3.6.**
   `lab/exp06_multiple_scattering.py` computes $a\,\delta n$ with the cytoplasm radius
   ($a=4.8\times0.62=2.98$ μm) and the *maximum* $|\delta n|$ in the volume (0.0769 at
   contrast 1). The volume-mean contrast of that phantom is 0.0215, a ratio of 0.28
   (computed, `realsample_table.py --phantom`). So the thresholds of
   [[multiple-scattering-check]] §5 and [[gridding-and-forward-model]] §4 are $0.23\lambda$ and $0.46\lambda$ in the
   max-contrast reading and $0.064\lambda$ and $0.13\lambda$ in the mean-contrast
   reading. [[chen1998]]'s $0.24\lambda$ is for uniform cylinders, where the two readings
   coincide. For a nearly uniform object (a red cell, a bacterium) the project's results
   therefore bracket the threshold between 0.06 and 0.24 and do not say where in that
   interval it lies *(inference)*.
2. **Single cells straddle the thresholds; nothing whole is safely below them.** Bacteria
   and isolated sub-micron organelles come out at 0.03–0.14. Red cells and yeast at
   0.1–0.5. Adherent mammalian cells at 0.07–0.45 by mean contrast and 0.15–0.7 by
   maximum contrast. Spheroids and tissue slices are above 1.
3. **GODT's printed phantom sits at about 0.15–0.30 through its height** (0.3–0.8 along
   its lateral semi-axes), using the indices and immersion oil of [[ziemczonok2019]].
   That is at or just above the point where the design margins have already shrunk to
   1.3× and 1.1×. The immersion liquid GODT itself used for the phantom was not found in
   the paper text, and the contrast depends entirely on it.
4. **GODT's neurons cannot be placed from the paper.** It gives no index, no cell
   thickness and no mounting medium that the retrieval found. Live cultured neurons in
   aqueous medium, taken as a proxy, give 0.6–1.2 by maximum contrast and 0.1–0.6 by an
   assumed mean contrast, which is at or beyond the upper threshold.
5. **So the accuracy margins of the design study describe a regime below most whole
   cells.** The geometric results (null comb, blind plane, the rule
   $\Delta\approx2\pi/K_{\max}$) do not depend on contrast
   ([[multiple-scattering-check]] §2) and apply to all rows.

---

## 1. Definitions

$x=a\,\delta n/\lambda$ with $\lambda$ the vacuum wavelength, $a$ a radius or half path
length and $\delta n$ the index contrast to the surrounding medium. Three choices are
hidden in it.

- **Which $\delta n$.** *max*: the largest contrast anywhere in the object, as exp06
  uses. *mean*: the volume average, which is what sets the phase accumulated through the
  object. Both are given where the source allows. Thresholds in the two readings, from
  [[gridding-and-forward-model]] §4 (multi-slice Born, checked against an exact solver;
  nearest-voxel / exact estimator) and the 0.28 ratio above:

  | exp06 contrast | $x_{\max}$ | $x_{\text{mean}}$ | fusion vs averaging | best pair vs single 0.8 μm |
  |---|---|---|---|---|
  | (linear model, no scattering) | n/a | n/a | 9.6× (on-voxel data); 1.8× / 5.9× (continuum data) | 14.3×; 2.5× / 8.8× |
  | 0.5 | 0.229 | 0.064 | 1.15× / 1.08× | 1.33× / 1.31× |
  | 1.0 | 0.458 | 0.128 | 1.02× / 0.96× | 1.09× / 1.07× |
  | 2.0 | 0.915 | 0.256 | ~1× (not validated) | ~1× (not validated) |

  No design margins were computed below contrast 0.5; exp06's linearity table has one
  lower point, $x_{\max}=0.114$ ($x_{\text{mean}}=0.032$), with 7 % nonlinearity.
  *(The thresholds were first measured with exp06's WPM-based model D,
  [[multiple-scattering-check]] §5, which gave the same margins at contrast 0.5.)*
- **Which $a$.** The phantom is a sphere, so radius and half path agree. Adherent cells
  and red cells are flat. The table gives the half thickness (normal-incidence path) and,
  where it matters, the lateral radius. At GODT's 60° illumination the path through a
  flat object is longer than its thickness, so the half-thickness figure is a lower
  bound *(inference)*.
- **Which $\lambda$.** GODT's experiment used a 530 nm LED and its simulation 0.5 μm
  ([[s42005-025-02416-3]]). The two columns differ by 6 %. Indices were measured at
  457–633 nm (one at 1310 nm) and are used unchanged. Dispersion of the *contrast* over
  that range is ignored *(inference; small against the ranges shown, except for the
  1310 nm row)*.

Medium indices used where a paper does not give its own: 1.333–1.338, spanning PBS
1.3332 and a buffer medium 1.3379 (Baczewska 2021) and 1.337 at 532 nm (Kim 2016).

## 2. The table

`uniform` means one reported index for the whole object, so max and mean coincide.
Ranges are input ranges multiplied through, not confidence intervals.

| sample | $a$ (μm) | $\delta n$ | medium | measured at | kind | $x$ at 530 nm | $x$ at 500 nm | source |
|---|---|---|---|---|---|---|---|---|
| **Reference: simulations** | | | | | | | | |
| exp06 phantom, contrast 1 | 2.98 | 0.0769 | 1.33 | n/a | max | 0.43 | 0.46 | computed |
| same | 2.98 | 0.0215 | 1.33 | n/a | mean | 0.12 | 0.13 | computed |
| GODT's Shepp–Logan, index 1.33–1.335, "maximum size" 40 μm | 20 (half the maximum size) | 0.005 | 1.33 | n/a | max | 0.19 | 0.20 | [[s42005-025-02416-3]] |
| **GODT's experimental samples** | | | | | | | | |
| Printed phantom, half height 5.5 μm, high-dose polymer 1.5167 | 5.5 | 0.029 | oil 1.4878 | 633 nm | max | 0.30 | 0.32 | [[ziemczonok2019]] |
| same, low-dose polymer 1.5021 | 5.5 | 0.014 | oil 1.4878 | 633 nm | mean, lower bound | 0.15 | 0.16 | [[ziemczonok2019]] |
| same, lateral semi-axes | 12.5–15 | 0.014–0.029 | oil 1.4878 | 633 nm | mean | 0.34–0.82 | 0.36–0.87 | [[ziemczonok2019]] |
| Printed "nucleolus", 2.6 μm, against the "nucleus" | 1.3 | 0.015 | polymer | 633 nm | uniform | 0.036 | 0.038 | [[ziemczonok2019]] |
| Neuron soma, RI 1.3424–1.3990 (proxy: live cortical culture) | 5–10 **assumed** | 0.061–0.066 | **assumed** 1.333–1.338 | 457 nm | max | 0.58–1.25 | 0.61–1.32 | Lee 2021 |
| same | 5–10 **assumed** | 0.01–0.03 **assumed** | same | n/a | mean | 0.09–0.57 | 0.10–0.60 | none; see §4 |
| **Other single cells** | | | | | | | | |
| HeLa, index 1.36–1.40 (Choi) or nucleolus 1.363 (Kim) | 3–6 **assumed** | 0.026–0.063 | **assumed** 1.337 | 633, 532 nm | max | 0.15–0.71 | 0.16–0.76 | Choi 2007; Kim 2019 |
| HeLa, cytoplasm 1.349 (Kim) to 1.36–1.39 (Choi) | 3–6 **assumed** | 0.012–0.038 | **assumed** 1.337 | 532, 633 nm | mean | 0.07–0.43 | 0.07–0.46 | Kim 2019; Choi 2007 |
| Red blood cell, 1.399 ± 0.006, half mean thickness | 0.93 | 0.066 | PBS, 1.3332 from Baczewska 2021 | 514 nm | uniform | 0.12 | 0.12 | Park 2008 |
| same, lateral radius | 4 | 0.066 | same | 514 nm | uniform | 0.50 | 0.53 | Park 2008 |
| Yeast, 1.411 (IQR 1.402–1.420), 25 fL (IQR 17–37) | 1.6–2.1 | 0.064–0.087 | **assumed** 1.333–1.338 | not retrieved | uniform | 0.19–0.34 | 0.20–0.36 | Midtvedt 2019 |
| *E. coli*, 1.384, relative index 1.0387, 0.5–1.4 fL | 0.49–0.69 | 0.052 | water | 589 nm | uniform | 0.05–0.07 | 0.05–0.07 | Balaev 2002; Oh 2020 |
| **Small, high-contrast structures** (contrast to their surroundings) | | | | | | | | |
| Lipid droplet, 1.386 against cytoplasm 1.347, 0.34–0.75 fL | 0.43–0.56 | 0.039 | cytoplasm | 532 nm | uniform | 0.03–0.04 | 0.03–0.04 | Kim 2016 |
| Lipid droplet, 1.468–1.479, size from Kim 2016, cytoplasm from Kim 2019 | 0.43–0.56 | 0.12–0.13 | cytoplasm | see source page | uniform | 0.10–0.14 | 0.10–0.15 | [[2024-06-14-598984-polarization-hilbert-holotomography]] |
| Nucleolus against nucleoplasm | 1–2 **assumed** | 0.01–0.03 | nucleoplasm | 532, 633 nm | uniform | 0.02–0.11 | 0.02–0.12 | Kim 2019; Choi 2007 |
| **Thick samples** | | | | | | | | |
| Liver-cell spheroid, diameter 140 μm | 70 | 0.01–0.03 **assumed** | not retrieved | not retrieved | mean | 1.3–4.0 | 1.4–4.2 | Yasuhiko & Takeuchi 2023 |
| Acute rat brain slice, grey matter 1.361–1.369, 300 μm | 150 | 0.019–0.027 | aCSF 1.342 ± 0.007 | 1310 nm | mean | 5–8 | 6–8 | Sun 2012 |

How the derived inputs were formed (all *(inference)* from the quoted source values):

- **Printed phantom.** Table 1 of [[ziemczonok2019]] gives the three polymer indices; the
  immersion oil is 1.4878 ± 0.0002. Which dose fills the bulk of the body was not stated
  in the passages read, so the body's mean contrast lies between the low-dose and
  high-dose rows. The nucleolus row is 1.5167 − 1.5021; the paper's design value is
  $\Delta$RI = 0.015.
- **Red cell.** Park 2008 gives volume 93.1 ± 7.9 fL and a projected diameter of about
  8 μm. Mean thickness is volume over disc area, 1.85 μm. Park 2008 names PBS without an
  index; 1.3332 is Baczewska 2021's PBS, whose wavelength was not found.
- **Yeast, *E. coli*, lipid droplets.** Radius of the sphere of the reported volume. An
  *E. coli* rod is longer and thinner than that sphere. The *E. coli* contrast is
  $1.384-1.384/1.0387$.
- **HeLa.** Choi 2007 reports tomogram slices at 5.5–12 μm above the substrate, a span of
  6.5 μm; the half thickness 3–6 μm is read from that. Its medium index was not in the
  document. Kim 2019 measured its medium with a refractometer but the value was not in
  the passages read. The mean row runs from Kim's cytoplasm value to the midpoint of
  Choi's cytoplasm range (1.375). The Choi and Kim cytoplasm values differ by 0.01–0.04,
  which is most of the width of the HeLa rows.
- **Nucleolus.** Kim 2019 counts only bodies above 2 μm as nucleoli, hence $a\ge1$ μm; the
  upper 2 μm is assumed.
- **Lipid droplets.** The two rows differ threefold in contrast. Kim 2016's 1.386 and the
  1.468–1.479 on the ingested source page are different cells, treatments and
  instruments, and this note does not reconcile them.

## 3. How each source was read

Three routes, in decreasing order of directness. **P**: passages of the full text
returned by a paper index and read directly. **Q**: the open full text was fetched by a
scraper and queried; verbatim sentences came back and only those were read. **A**:
abstract only.

| source | route | what was taken |
|---|---|---|
| [[s42005-025-02416-3]] (ingested) | Q | Shepp–Logan index range and size, 530 nm LED, NA 1.45, neuron culture protocol. No phantom immersion, neuron mounting medium or index values were returned |
| [[ziemczonok2019]] (ingested) | P, Q | dimensions, Table 1 indices (P); oil index, 632.8 nm (Q) |
| Choi 2007 | Q (author-hosted PDF) | HeLa compartment indices, slice heights, 633 nm |
| Kim 2019 | P | HeLa nucleolus 1.363, cytoplasm and nucleoplasm 1.349, 532 nm |
| Kim 2016 | P, Q | medium 1.337 at 532 nm (P); droplet and cytoplasm indices, droplet volume (Q) |
| Park 2008 | Q | red cell index, volume, diameter, 514 nm |
| Midtvedt 2019 | P | yeast index and volume |
| Oh 2020 | P | *E. coli* volume range; no index value in the passages |
| Balaev 2002 | A | *E. coli* index, from a model fit to suspension extinction data, not a single-cell measurement |
| Lee 2021 | Q (group-hosted PDF) | soma and neurite index ranges, 457 nm |
| Baczewska 2021 | Q | medium and PBS indices; fixation lowers RI, most in the nucleolus. No numeric cell indices were returned (they are in a figure) |
| Yasuhiko & Takeuchi 2023 | A, P | spheroid diameter; segmentation thresholds $n>1.340$ (spheroid) and $n>1.375$ (lipid) |
| Sun 2012 | Q | grey-matter and aCSF indices, slice thickness, 1310 nm OCT |

Searched and not usable: Rappaz et al. 2005 (neuron index and thickness by digital
holography; only the abstract was reachable and it has no numbers) and Liu et al. 2014
(single-bacterium refractometry; abstract without values).

## 4. Reading the table against the thresholds

**Well below the tested range: bacteria and isolated small organelles.** *E. coli* at
0.05–0.07 and the Kim 2016 lipid droplet at 0.03–0.04 are below exp06's lowest contrast
in the max reading. This is where the linear-model margins would apply with the least
reduction. Two limits. exp06 did not run the design comparison there, so "undiminished"
is an extrapolation of its trend, not a result *(inference)*. And an organelle is not
imaged alone: it sits inside a cell whose own $x$ sets how far the incident field is
already distorted when it reaches the organelle. The organelle rows apply to an isolated
droplet or to one in a thin, weak cell region.

**Around the first threshold: red cells through their thickness, yeast, thin cells by
mean contrast.** 0.1–0.35. For these nearly uniform objects headline 1 applies: the
threshold is somewhere in 0.06–0.24, and they are inside or just above that interval.
The margins to expect are of the order of the 1.3× and 1.1× measured at the first
threshold, not the 14.3× and 9.6× of the linear model.

**At or beyond the second threshold: adherent mammalian cells by maximum contrast, red
cells along their radius.** 0.15–0.7 for HeLa and, by proxy, 0.6–1.2 for neuron somas.
[[gridding-and-forward-model]] §4 finds the design margin gone at 0.46 with a validated
forward model, and exp06 found the first-order estimate no better than zero by 0.92
([[multiple-scattering-check]] §5). A nonlinear inverse of the kind [[ma2017]]
describes is what both notes name as the fix.

**Far beyond: spheroids and tissue slices.** $x>1$ under any reasonable contrast. The
spheroid paper is itself a method for getting past the linear model in such samples.

**GODT's two samples.**

- *Printed phantom.* 0.30 (max) and at least 0.15 (mean) through the height, both above
  the first threshold in their own reading (0.23 and 0.064) and at or above the second
  in the mean reading (0.13). On this evidence the phantom is in the regime where exp06
  finds margins of 1.0–1.3×. This holds for the oil of [[ziemczonok2019]] and indices at
  633 nm; GODT imaged it at 530 nm in a liquid this note could not identify.
- *Hippocampal neurons.* Not determinable. The cells are fixed, and Baczewska 2021
  reports that fixation lowers cell RI, so live-cell values are an upper proxy at best.
  A mounting medium other than buffer would change the contrast by more than everything
  else in the row. If the sample resembles live cultured neurons in aqueous medium, the
  soma is at or beyond the second threshold; neurites, at 1.3285–1.3424 in Lee 2021 and
  about a micron across *(size assumed)*, would be far below the first.

**What this means for the design study.** Its accuracy margins are properties of weakly
scattering objects, and most whole cells are not in that class by this axis. Three
things soften that. The geometric results hold at any contrast. Most of the reduction
from 9.6× to 1.8× is gridding, which is fixable and independent of the sample: without
it the margin on single-scattering data is 5.9× ([[gridding-and-forward-model]]). And
the multiple-scattering margins rest on one phantom and a forward model validated only
on near-axis rays and up to $x_{\max}=0.46$ ([[full-wave-referee]] §3). An inverse that
models multiple scattering would move the thresholds; none has been tested.

## 5. What this does *not* establish

- **No sample was measured or simulated.** The $x$ column is arithmetic on literature
  values; the thresholds it is compared with come from one phantom at one size, with
  $x$ varied through contrast only ([[multiple-scattering-check]] §7).
- **Where a uniform object's threshold lies.** The factor 3.6 between the two readings is
  unresolved for uniform objects. A uniform-sphere run of exp06 would settle it.
- **GODT's actual operating point.** Phantom immersion, neuron mounting medium and neuron
  thickness were not found. The published paper's Methods were queried, not read end to
  end, and its supplementary material was not checked.
- **Sizes marked assumed** (neuron soma radius, HeLa half thickness, nucleolus radius,
  spheroid contrast, neuron mean contrast) have no retrieved source.
- **Route Q is one step removed.** Sentences came back as quotes from a query over the
  page. A sentence the query missed would not have been seen, and "not returned" is
  weaker than "not in the paper".
- **Literature RI values are themselves first-order reconstructions.** Most come from
  Rytov or projection tomography with a missing cone, which underestimates compact
  high-index features ([[ziemczonok2019]] reports an 8 % cumulative RI loss on its
  spheres). The lipid-droplet disagreement may be an instance. Maximum-contrast rows are
  more likely low than high *(inference)*.
- **Not a survey.** One or two papers per sample type, chosen for reachable full text.
  Cell-to-cell and cell-cycle variation is represented only where a paper gave a spread.
- **$a\,\delta n/\lambda$ is one number for a 3D object.** It ignores shape, internal
  structure and the 60° illumination path, each of which changes the scattering regime
  at fixed $x$.

## Reproduction

```bash
cd lab
python realsample_table.py             # the x column, from the inputs in the file
python realsample_table.py --phantom   # exp06 phantom: dn_max, mean, ratio 0.28
```

## References (web-sourced, not ingested)

- Balaev, Dvoretski & Doubrovski, "Refractive index of escherichia coli cells," *Proc. SPIE* **4707**, 253–260 (2002). doi:10.1117/12.475627
- Baczewska, Eder, Ketelhut, Kemper & Kujawińska, "Refractive index changes of cells and cellular compartments upon paraformaldehyde fixation acquired by tomographic phase microscopy," *Cytometry A* **99**, 388–398 (2021). doi:10.1002/cyto.a.24229. Ref. 64 of [[s42005-025-02416-3]].
- Choi et al., "Tomographic phase microscopy," *Nat. Methods* **4**, 717–719 (2007). doi:10.1038/nmeth1078
- T. K. Kim et al., "Physicochemical properties of nucleoli in live cells analyzed by label-free optical diffraction tomography," *Cells* **8**(7), 699 (2019). doi:10.3390/cells8070699
- K. Kim, Lee, Yoon, Heo, Choi & Park, "Three-dimensional label-free imaging and quantification of lipid droplets in live hepatocytes," *Sci. Rep.* **6**, 36815 (2016). doi:10.1038/srep36815
- Lee, Hugonnet, Nam, Yoon, Park, Park, Han & Park, "Label-free monitoring of 3D cortical neuronal growth in vitro using optical diffraction tomography," *Biomed. Opt. Express* **12**(11), 6928–6944 (2021). doi:10.1364/BOE.439817. Author order as returned by the PDF query.
- Midtvedt, Olsén, Höök & Jeffries, "Label-free spatio-temporal monitoring of cytosolic mass, osmolarity, and volume in living cells," *Nat. Commun.* (2019). doi:10.1038/s41467-018-08207-5
- Oh et al., "Three-dimensional label-free observation of individual bacteria upon antibiotic treatment using optical diffraction tomography," *Biomed. Opt. Express* **11**, 1257–1267 (2020). doi:10.1364/BOE.377740. Ref. 65 of [[s42005-025-02416-3]].
- Park, Diez-Silva, Popescu, Lykotrafitis, Choi, Feld & Suresh, "Refractive index maps and membrane dynamics of human red blood cells parasitized by *Plasmodium falciparum*," *PNAS* **105**(37), 13730 (2008). doi:10.1073/pnas.0806100105
- Sun, Lee, Wu, Sarntinoranont & Xie, "Refractive index measurement of acute rat brain tissue slices using optical coherence tomography," *Opt. Express* **20**(2), 1084–1095 (2012). doi:10.1364/OE.20.001084
- Yasuhiko & Takeuchi, "In-silico clearing approach for deep refractive index tomography by partial reconstruction and wave-backpropagation" (2023). doi:10.1038/s41377-023-01144-z
- Not usable: Rappaz et al., *Opt. Express* **13**(23), 9361 (2005), doi:10.1364/opex.13.009361; Liu et al., *Lab Chip* **14**(21) (2014), doi:10.1039/c4lc00783b.
