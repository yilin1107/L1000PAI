# Chem family taxonomy v2 curation policy

Version: 2.0.0  
Lock date: 2026-08-26

## Scope

The v2 term registry covers every exact compound term in the paired LINCS Chem
consensus library. Pharmacologic family membership is curated only for the 80
terms used by the historical ten-case Chem-family analysis after adding the
omitted `PLX-4720` input signature. All other terms remain explicitly
`unreviewed`; absence of a family membership is not evidence that a compound
has no mechanism.

## Evidence boundary

The historical ten-case roster was assembled after its source result tables
existed. It is therefore labelled `retrospective_author_curated`, not
`prespecified`. Taxonomy v2 is frozen for future analyses from the lock date;
legacy outcome columns were not copied into any frozen taxonomy table.

## Family levels

- `direct_moa`: compounds share a direct pharmacologic mechanism and may enter
  a strict same-mechanism analysis when the case mapping selects that family.
- `pathway_context`: compounds perturb a common signaling or phenotypic axis,
  but do not necessarily share the same direct target. This level is
  exploratory and must not be described as same-mechanism evidence.
- `phenotypic_class`: reserved for future compound classes defined primarily by
  a shared phenotype rather than direct target pharmacology.

## Identity and aggregation

- `library_term` is the exact primary key. Punctuation-normalized text is only
  a candidate-matching aid and is never used to merge compounds automatically.
- Every library term maps to a stable `analysis_unit_id`. Family denominators
  use one vote per analysis unit.
- `17-AAG` and `Tanespimycin` share one chemical identity and analysis unit.
- `Ixazomib` and `MLN-9708` retain distinct chemical-form identifiers but share
  one drug-program analysis unit.
- Raw term-level scores and adjusted p-values must remain available. A display
  or family summary may select the frozen preferred term, but must not choose a
  representative term based on the observed activity score.

## Direction semantics

For Chem paired signatures, a positive L1000PAI score denotes transcriptional
concordance with the compound-induced signature and a negative score denotes
opposition. Family annotation does not by itself establish drug sensitivity,
target engagement, efficacy, or causal dependency.

## Change control

Any change to a family definition, membership, alias relationship, or preferred
analysis term requires a new taxonomy version and regenerated SHA-256 manifest.
Historical versions remain immutable.
