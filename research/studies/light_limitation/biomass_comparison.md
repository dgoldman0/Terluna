# Annual biomass production compared with Earth

31 October 2026. A provisional comparison for research planning, using the ecology branch's
production products and the current climate programme's status. The reference year is 365.25
Earth days. Production means net primary production (NPP): carbon incorporated into new biomass
after the producers' own respiration. It supplies food webs and tissue replacement; annual
production can greatly exceed the change in standing biomass or long-term carbon storage.

**About 20 Gt C a year, with a working range of 15–25, is the present estimate for the Open
Moon's sunlight-supported biosphere.** The range is a rounded judgment while the revised 301 K
climate and biological responses remain unresolved. It is about 14–24% of the historical Earth
benchmark of 104.9 Gt C a year, or 13–25% using a rounded Earth comparison of 100–120. A world
with 7.44% of Earth's surface would therefore average roughly two to three times Earth's
production per square metre across land and water together. Assuming dry biomass is 45% carbon
by mass in both worlds, the lunar estimate corresponds to 33–56 billion tonnes of dry biomass
a year, with about 44 billion at the working centre. These are conditional capabilities of an
established engineered biosphere.

The numerical basis is [biomass_comparison.py](biomass_comparison.py), its explicit
[inputs and source records](biomass_comparison_inputs.json), and the generated
[comparison](results/biomass_comparison.json). The calculation retains the existing land
production, reconciles the waters with the nutrient-poor design, and reports the sky separately.
Its low and high cases are choices of parameters, not statistical limits.

| Component | Annual production, Gt C | Evidence and conditions |
|---|---:|---|
| Rain-fed dry land | 14.0–14.6 | Existing screen on the older weather; engineered canopy carbon budget with water limits |
| Terrestrial wet margins | 1.17–5.44 | Existing 0.49–2.01 million km² scenario with ample water and nutrients |
| Nutrient-poor lakes and open seas, plus productive coasts | 0.26–2.44 | Reconciled scenario described below |
| Surface total | 15.4–22.5 | Land and waters; illustrative middle 18.8 |
| Aerophyte host production | 0.02–0.83 | At 0.1–2% cover; shown separately until displaced surface production is known |

The dry land is 23.4 million km². The food study already removes both the common-level seas
(28% of the Moon, 10.62 million km²) and the additional rain-fed lakes (3.92 million km²).
Terrestrial wet margins replace rain-fed land within that allocation; they are not extra area.
Sky organisms can use light above clouds or above otherwise unproductive ground, but their
interception also reduces the light below them. Their host production therefore cannot simply
be added to an unchanged surface biosphere. The rounded 15–25 planning range allows discussion
of these unresolved effects without declaring its endpoints validated bounds.

The older [light-limitation calculation](README.md) gave 16.2–25.5 Gt C/year. Its lake rates,
100–400 g C/m²/year, represented mesotrophic to eutrophic conditions. The later aquatic
direction keeps open lakes and seas nutrient-poor, concentrating richer communities in margins,
shelves, deltas and farms. This comparison uses the oxygen study's oligotrophic gross production
bracket of 50–300 mg C/m²/day, central 200, as recorded secondarily from Wetzel. At 200 mg,
annual GPP is 73.05 g C/m². The oxygen model prescribes this whole-cycle mean and calibrates
its light response to achieve it; the rate is an assumption about lunar organisms. An explicitly assumed NPP/GPP fraction of 0.4–0.6 converts gross
production into new biomass; its middle gives 36.5 g C/m²/year. The oxygen study's euphotic
community respiration includes consumers and decomposers, so it cannot supply that fraction.

The constructed aquatic cases apply the same oligotrophic rates to all lake area and to the
open part of the seas, retaining productive coasts at the old screen's guessed 300–1,000
g C/m²/year over 5–15% of sea area. This gives 0.26–2.44 Gt C/year, with 1.18 in the illustrative
middle. A separately allocated productive aquatic lake margin remains absent. Coast and open-sea
fractions sum to one in each case. The earlier screen paired their endpoints independently,
giving combined sea areas of 90% and 110% of the actual area; its historical result remains
unchanged, while this comparison uses complementary areas.

The completed October 20 climate case, `A28R_301k_s0897`, reached 301.25 K at sunlight scale
0.897. Subsequent changes to convection, vertical resolution, surface exchange, clouds and
precipitation alter that tuning. The October 30 synthesis anticipated approximately 3.5–4.5%
more sunlight before the final rain and snow treatment was tuned, implying scale 0.928–0.937.
The revised configuration still needs its short tuning and long design runs. Its snow smoke
test establishes implementation behaviour, not equilibrium weather.

The biomass calculation's weather is `A28_dim5_moon`; both that weather and its source spectra
use scale 0.950. The later old climate lock
at 0.970 is a different reference. Relative to 0.950, the provisional revised sunlight is only
1.3–2.3% lower, versus 5.6% lower in the completed 0.897 case. If surface PAR changed in the same
proportion and everything else stayed fixed, the old screen's whole-biosphere light elasticity
of 0.07–0.295, whose weighting retains the older water and sky mixture, would imply a direct NPP reduction of about 0.1–0.7% in the revised case, versus
0.4–1.6% in the completed case. These are sensitivities, not corrections applied to the totals:
the clouds, spectral transfer, temperatures and rain have also changed, and the final sunlight
scale remains open. Climate source identities are recorded in the inputs file, including the
read-only main-branch synthesis at commit `2f6cca9`.

Water and respiration matter more to this estimate than that small direct light change. The
older interpolation's 25% increase in mean land rain at 301 K is insufficient evidence for a
25% biomass increase: rain distribution, soil water, transpiration demand, warmer respiration
and acclimation all affect the result. The pending climate should supply the temperature and
water fields for the next canopy and land calculation. Nitrogen and phosphorus sufficiency,
continuous-light tolerance, long-dark survival and complete ecological life cycles remain
conditions of the estimate. The plant budgets explicitly include the long lunar night, while aquatic daily rates are
prescribed averages over the whole cycle. Neither warrants another factor of one half; the
aquatic assumption still needs a physiological test of night survival and daytime production.

The large amount of land, warm growing conditions and efficient use of diffuse light support
substantial production per unit planetary area in the existing models. Low gravity allows
different architectures but supplies no additional photosynthetic energy. With water available
everywhere, the older land calculation gives 46.8–54.5 Gt C/year; including the aquatic cases
gives approximately 47–57 Gt C/year, around half Earth's total. That is a conditional resource
potential on the older climate, requiring ample nutrients, adapted organisms and water delivery.
The changed vegetation's albedo and water cycle would themselves require another climate
calculation. Powered microbial food production and electrically lit cultivation are outside
the sunlight-supported totals here, and none of the totals specifies a sustainable harvest.

For Earth, [Field et al. (1998)](https://pubmed.ncbi.nlm.nih.gov/9657713/) report 104.9 Pg C/year;
[Behrenfeld et al. (2001)](https://pubmed.ncbi.nlm.nih.gov/11283369/) report 111–117 during
1997–2000. They support the rounded comparison scale, rather than a claim to know Earth's 2026
total. The more recent [Zhang et al. (2025)](https://www.nature.com/articles/s41558-025-02375-1)
find different production trends on land and in the ocean. Abstracts were checked for this note;
no full scholarly source review or empirical validation of lunar organisms is claimed.

Reproduce the arithmetic with `python -m research.studies.light_limitation.biomass_comparison`;
run its checks with `python -m pytest research/studies/light_limitation/test_biomass_comparison.py`.
The checks exercise area conservation, conversion of carbon units, the GPP/NPP distinction,
separate sky reporting and regeneration of the recorded result. They do not validate the climate
or biological assumptions.

Verification on 31 October 2026: `make check` passed, including 2,825 Python tests
with 20 skips, the JavaScript and calendar checks, provenance checks and workspace
validator. The independent comparison review found no blocking arithmetic or provenance
issues. These checks establish repository consistency, not biological validation.
