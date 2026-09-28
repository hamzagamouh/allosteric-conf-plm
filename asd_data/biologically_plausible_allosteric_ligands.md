# Biologically plausible allosteric ligands in `asd_data/*/allosteric_ligands_*.json`

**Scope.** The 10 fold files (`train/` + `val/`) describe 270 unique proteins.
Each protein carries an `"Allosteric ligands"` list whose entries have the form
`pdb-chain-CCD-resnum`. Extracting the 3‑letter PDB chemical‑component code (CCD)
from every entry gives **1 697 distinct ligand codes** bound at annotated
allosteric sites across all related PDB structures.

Names were resolved from the RCSB Chemical Component Dictionary
(`data.rcsb.org`). The overwhelming majority of the 1 697 codes are synthetic
medicinal‑chemistry compounds (long IUPAC names, drug code names such as
*Cobimetinib*, *Trametinib*, *Vandetanib*, *Selpercatinib*, *Rapamycin*
analogues, kinase‑inhibitor series `X…`, `Y…`, `1J…`, `7U…`, etc.). Those are
**excluded** below.

What remains — molecules that plausibly occur in a biological system
(endogenous metabolites, cofactors, ions, hormones, neurotransmitters,
signalling molecules, lipids, sugars, and microbial/plant natural products) —
is listed here. The number in parentheses is the count of the 270 proteins in
which that code appears at an allosteric site.

---

## 1. Nucleotides, nucleosides, nucleobases (endogenous)

| Code | Name | (proteins) |
|------|------|-----------|
| ATP | adenosine‑5′‑triphosphate | 16 |
| ADP | adenosine‑5′‑diphosphate | 18 |
| AMP | adenosine‑5′‑monophosphate | 6 |
| CMP | adenosine‑3′,5′‑cyclic monophosphate (cAMP) | 6 |
| ACK | 2′,3′‑cyclic AMP | 1 |
| ADN | adenosine | 3 |
| ADE | adenine | 2 |
| MTA | 5′‑deoxy‑5′‑methylthioadenosine | 1 |
| AR6 | ADP‑ribose | 1 |
| GTP | guanosine‑5′‑triphosphate | 5 |
| 5GP | guanosine‑5′‑monophosphate (GMP) | 2 |
| PCG | cyclic GMP (cGMP) | 4 |
| 35G | guanosine‑3′,5′‑monophosphate | 1 |
| GUN | guanine | 1 |
| HPA | hypoxanthine | 1 |
| CTP | cytidine‑5′‑triphosphate | 3 |
| TTP | thymidine‑5′‑triphosphate | 3 |
| TMP | thymidine‑5′‑monophosphate | 1 |
| UDP | uridine‑5′‑diphosphate | 1 |
| DTP | 2′‑deoxyadenosine‑5′‑triphosphate (dATP) | 4 |
| DGT | 2′‑deoxyguanosine‑5′‑triphosphate (dGTP) | 3 |
| DCP | 2′‑deoxycytidine‑5′‑triphosphate (dCTP) | 1 |
| DCM | 2′‑deoxycytidine‑5′‑monophosphate (dCMP) | 1 |
| DU  | 2′‑deoxyuridine‑5′‑monophosphate (dUMP) | 1 |
| DUT | 2′‑deoxyuridine‑5′‑triphosphate (dUTP) | 1 |
| U2G | uridylyl‑2′,5′‑guanosine (dinucleotide) | 1 |

## 2. Second messengers / nucleotide alarmones

| Code | Name | (proteins) |
|------|------|-----------|
| C2E | cyclic di‑GMP (c‑di‑GMP) | 2 |
| 2BA | cyclic di‑AMP (c‑di‑AMP) | 1 |
| G4P | guanosine‑5′,3′‑tetraphosphate (ppGpp) | 1 |
| 0O2 | guanosine‑5′‑triphosphate‑3′‑diphosphate (pppGpp) | 1 |
| I3P | D‑myo‑inositol‑1,4,5‑trisphosphate (IP₃) | 1 |
| 4IP | inositol‑1,3,4,5‑tetrakisphosphate (IP₄) | 2 |
| IHP | inositol hexakisphosphate (InsP₆ / phytic acid) | 3 |
| NO  | nitric oxide | 1 |
| PEO | hydrogen peroxide | 1 |

## 3. Redox / acyl / methyl cofactors (endogenous)

| Code | Name | (proteins) |
|------|------|-----------|
| NAD | NAD⁺ | 3 |
| NAI | NADH | 2 |
| NAP | NADP⁺ | 1 |
| NDP | NADPH | 2 |
| NMN | nicotinamide mononucleotide (NMN) | 1 |
| NCN | nicotinate mononucleotide | 1 |
| SAM | S‑adenosyl‑L‑methionine | 2 |
| SAH | S‑adenosyl‑L‑homocysteine | 2 |
| COA | coenzyme A | 1 |
| ACO | acetyl‑coenzyme A | 1 |
| PKZ | palmitoyl‑CoA | 1 |
| CMC | carboxymethyl‑CoA | 1 |
| PLP | pyridoxal‑5′‑phosphate (vitamin B6 cofactor) | 1 |
| BTN | biotin (vitamin B7) | 1 |
| BTX | biotinyl‑5′‑AMP | 1 |
| PAU | pantothenic acid (vitamin B5) | 1 |
| PAZ | 4′‑phosphopantothenate | 1 |
| FOL | folic acid | 1 |
| DHF | dihydrofolate | 1 |
| THG | (6S)‑5,6,7,8‑tetrahydrofolate | 1 |
| FFO | folinic acid (5‑formyl‑THF) | 1 |
| HBI | 7,8‑dihydrobiopterin | 1 |

## 4. Sugar‑phosphates & central‑carbon metabolites

| Code | Name | (proteins) |
|------|------|-----------|
| G6P | glucose‑6‑phosphate (α) | 2 |
| BG6 | glucose‑6‑phosphate (β) | 2 |
| G16 | glucose‑1,6‑bisphosphate | 1 |
| FBP | fructose‑1,6‑bisphosphate | 1 |
| G3P | sn‑glycerol‑3‑phosphate | 2 |
| G3H | glyceraldehyde‑3‑phosphate | 1 |
| 13P | dihydroxyacetone phosphate | 1 |
| R5P | ribose‑5‑phosphate | 1 |
| M6D | mannose‑6‑phosphate (β) | 1 |
| GLP | glucosamine‑6‑phosphate (α) | 1 |
| 4R1 | glucosamine‑6‑phosphate (β) | 1 |
| 16G | N‑acetylglucosamine‑6‑phosphate (α) | 2 |
| 4QY | N‑acetylglucosamine‑6‑phosphate (β) | 1 |
| 0NZ | 2‑deoxyglucose‑6‑phosphate | 1 |
| PYR | pyruvate | 3 |
| 2OP | L‑lactate | 1 |
| PPY | phenylpyruvate | 1 |
| KPV | 5‑phenyl‑2‑oxovalerate | 1 |
| CIT / FLC | citrate | 3 / 1 |
| ICT | isocitrate | 1 |
| AKG | 2‑oxoglutarate (α‑ketoglutarate) | 3 |
| MLI / MLA | malonate / malonic acid | 3 / 1 |
| GOA | glycolate | 1 |
| PRE | prephenate | 1 |
| IGP | indole‑3‑glycerol phosphate | 1 |
| 1AL | allantoate (purine catabolism) | 1 |
| SHF | levulinate (5‑oxopentanoate) | 1 |
| POP | inorganic pyrophosphate (PPi) | 1 |

## 5. Sugar–nucleotides & glycans

| Code | Name | (proteins) |
|------|------|-----------|
| UD1 | UDP‑N‑acetylglucosamine | 2 |
| UDX | UDP‑xylose | 1 |
| GDD | GDP‑α‑D‑mannose | 1 |
| NCC | CMP‑N‑acetylneuraminic acid (CMP‑sialic acid) | 1 |
| BGC / GLC | β‑ / α‑D‑glucose | 2 / 3 |
| MAN | α‑D‑mannose | 1 |
| NAG | N‑acetyl‑β‑D‑glucosamine | 2 |
| BM7 | N‑acetyl‑β‑D‑mannosamine | 1 |
| INS | myo‑inositol | 1 |

## 6. Amino acids & physiological derivatives

Standard amino acids seen as allosteric ligands:
**ARG (5), TYR (5), PHE (4), TRP (6), GLN (3), GLU (2), HIS (2), SER (2),
CYS (2), ASP (1), GLY (1), ILE (1), LEU (1), LYS (1), MET (1), PRO (1),
THR (1)**, plus **DTY** (D‑tyrosine, 1).

Physiological modifications / related:

| Code | Name | (proteins) |
|------|------|-----------|
| SEP | phosphoserine | 4 |
| TPO | phosphothreonine | 2 |
| PTR | phosphotyrosine | 2 |
| MLY / MLZ / M3L | di‑ / mono‑ / tri‑methyl‑lysine | 3 / 1 / 1 |
| PCA | pyroglutamate | 1 |
| OAS | O‑acetyl‑L‑serine | 1 |
| NLG | N‑acetyl‑L‑glutamate | 1 |
| GGL | γ‑L‑glutamate | 1 |
| HRG | L‑homoarginine | 1 |
| ABU | γ‑aminobutyric acid (GABA) | 1 |

## 7. Amines, polyamines, neurotransmitters

| Code | Name | (proteins) |
|------|------|-----------|
| PUT | putrescine (1,4‑diaminobutane) | 2 |
| SPD | spermidine | 2 |
| SPM | spermine | 2 |
| N2P | cadaverine (1,5‑diaminopentane) | 1 |
| 13D | 1,3‑diaminopropane | 1 |
| SP5 | N¹‑acetylspermine | 1 |
| HLG | N¹‑acetylspermidine | 1 |
| SRO | serotonin | 1 |
| ACH | acetylcholine | 2 |
| GAI | guanidine | 2 |
| URE | urea | 2 |

## 8. Fatty acids & lipids

| Code | Name | (proteins) |
|------|------|-----------|
| PLM | palmitic acid (C16:0) | 4 |
| MYR | myristic acid (C14:0) | 4 |
| STE | stearic acid (C18:0) | 1 |
| OLA | oleic acid (C18:1) | 5 |
| HXA | docosahexaenoic acid (DHA, C22:6) | 3 |
| OCA | octanoic acid (C8:0) | 1 |
| 6NA | hexanoic acid (C6:0) | 1 |
| F15 | pentadecanoic acid (C15:0) | 1 |
| FPP | farnesyl diphosphate | 1 |
| PTG | 15‑deoxy‑Δ¹²,¹⁴‑prostaglandin J₂ | 1 |
| CLR | cholesterol | 7 |
| 5JK | 7α‑hydroxycholesterol (oxysterol) | 1 |
| Membrane phospholipids (mostly short‑chain synthetic surrogates of endogenous lipids): PC1, PC8, PX4, MC3, PLC, LPC, LPX, PEE, PEX, PGW, POV, PTY, PIO, DDR, NKP, OLC | – | – |

## 9. Hormones & signalling steroids/eicosanoids

| Code | Name | (proteins) |
|------|------|-----------|
| HCY | cortisol (hydrocortisone) | 1 |
| C0R | corticosterone | 1 |
| STR | progesterone | 1 |
| TES | testosterone | 1 |
| AS4 | aldosterone | 1 |
| T44 | 3,3′,5,5′‑tetraiodo‑L‑thyronine (thyroxine, T4) | 1 |
| T4A | 3,3′,5,5′‑tetraiodothyroacetic acid (TRIAC, T4 metabolite) | 1 |
| REA | all‑trans retinoic acid | 1 |
| 9CR | 9‑cis‑retinoic acid | 1 |
| GA4 | gibberellin A4 (plant hormone) | 1 |

*(Dozens of additional synthetic retinoid / rexinoid analogues — 2VP, 2VZ, 2W0,
2VR, 3RB, 3SW, 3T2, 754, 3TN, 4TN, 5TN, L79, BM6, LG2, LG3, 29V, R4M, OI2, OI5,
2QO, 2E3, 9HF, 9RA, 4XW … — are **excluded** as synthesized compounds.)*

## 10. Inorganic ions & small inorganic species (biologically essential)

| Code | Name | (proteins) |
|------|------|-----------|
| MG | Mg²⁺ | 51 |
| CA | Ca²⁺ | 28 |
| ZN | Zn²⁺ | 13 |
| MN | Mn²⁺ | 8 |
| FE | Fe³⁺ | 2 |
| CU / CU1 | Cu²⁺ / Cu⁺ | 2 / 4 |
| NA | Na⁺ | 26 |
| K | K⁺ | 8 |
| CO | Co²⁺ | 2 |
| NI | Ni²⁺ | 3 |
| PO4 | phosphate | 27 |
| SO4 | sulfate | 66 |
| SO3 | sulfite | 1 |
| BCT | bicarbonate | 2 |
| NO3 | nitrate | 2 |
| NO2 | nitrite | 1 |
| NH4 | ammonium | 1 |

*(SO4, PO4, NA, CL are also the commonest crystallization additives — see notes.)*

## 11. Microbial / plant natural products (not de‑novo synthetic drugs)

| Code | Name | (proteins) |
|------|------|-----------|
| STU | staurosporine | 2 |
| RAP | rapamycin (sirolimus) | 1 |
| FK5 | ascomycin (FK506 analogue) | 1 |
| FUA | fusidic acid | 1 |
| GET | geneticin / G418 (aminoglycoside) | 1 |
| IVM | ivermectin | 2 |
| SFG | sinefungin | 1 |
| GNT | (−)‑galantamine | 1 |
| HUP | huperzine A | 1 |
| HUB | huperzine B | 1 |
| NCT | nicotine | 1 |
| EPJ | epibatidine | 1 |
| L0B | lobeline | 1 |
| P0T | cannabidiol | 2 |
| AFT | aflatoxin B1 | 1 |
| PBQ | pentabromopseudilin | 1 |
| SJA | amorphadiene | 1 |
| S1A | soraphen A | 1 |
| LQ4 | piperlongumine | 1 |
| ROA | rosmarinic acid | 1 |
| FER | ferulic acid | 1 |
| HC4 | p‑coumaric acid (4‑hydroxycinnamic acid) | 1 |

---

## Notes / judgement calls

* **Non‑hydrolysable analogues of endogenous nucleotides** were kept separate
  and are *not* in the lists above, though they mimic natural ligands:
  ANP (AMP‑PNP, 12), AGS (ATP‑γ‑S, 2), ACP (AMP‑PCP, 2), APC (AMP‑PCP, 2),
  GN3 (GppNHp, 1), 75G / SP1 / RP1 / RP2 (cyclic‑nucleotide phosphorothioates),
  AN2 (AMP‑phosphoramidate), 2PN (imidodiphosphate), 3PO (triphosphate),
  DG3 (ddGTP), PMP (amino‑PLP), HF7/HDV/HEJ/HFD/6AT/HF4/0KX/F6G/GTF (nucleotide
  drug/analogue triphosphates).
* **Phosphate‑mimic ions** BEF (BeF₃⁻, 3), AF3 (AlF₄⁻, 1), MGF (MgF₃⁻, 1) are
  transition‑state analogues, not natural ligands.
* **Common crystallization / buffer additives** are not "synthesized drugs" but
  are not biologically meaningful modulators either: GOL/EDO/PEG/PGE/PG4/P6G/
  1PE/2PE/PE4/P3G/1PG/AE4/FWN (glycols & PEGs), MPD/MRD, DMS/DMSO, ACT/ACY/FMT
  (acetate/formate), EOH/IPA/PEL (alcohols), MES/EPE/BTB/TRS/TAM/CXS/CAC/MPD
  (buffers), IOD/BR/CS/CD/CL/AU/AG/HG/XE (heavy‑atom / non‑physiological ions),
  DOD (D₂O), MSE (selenomethionine – labelling), LMT/MA4/BOG/LMD/LSM/DDQ/CPS
  (detergents), UNX/UNL/UNK (unknown density). Treat CL, NA, SO4, PO4 as
  dual‑use (both physiological and additive).
* Amino‑acid codes also appear as ordinary residues; here they are counted only
  where listed as an `"Allosteric ligands"` entry (i.e. a free amino‑acid
  effector, e.g. Phe/Tyr/Trp feedback inhibition of DAHP synthase, Arg/Gln
  regulation, etc.).
