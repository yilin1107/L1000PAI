# Licence of the Chem taxonomy 3.0.1 data files

The files in this folder are data, not code. The MIT licence of the L1000PAI
package code does not apply to them.

They are released under the Creative Commons Attribution-ShareAlike 4.0
International licence (CC BY-SA 4.0), because they are adapted from sources
under share-alike terms:

| Source | Used for | Terms |
|---|---|---|
| ChEMBL (release 37) | molecule hierarchy, drug mechanisms, bioactivities, target classification | CC BY-SA 3.0 |
| IUPHAR/BPS Guide to PHARMACOLOGY (2026.3.1) | primary targets and affinities | database ODbL; contents CC BY-SA 4.0 |
| HGNC | gene symbols and UniProt accessions | CC0 |
| LINCS 2020 compound metadata (CLUE, Broad Institute) | compound identifiers, structure keys and the derived target/action call | public LINCS2020 input verified; see source-specific policy and attribution below; no source text is included |

When you use these files, please cite ChEMBL (Zdrazil et al., Nucleic Acids
Res 2024) and the Guide to PHARMACOLOGY (Harding et al., Nucleic Acids Res
2024), as well as L1000PAI.

`chem_taxonomy_v3_manifest.csv` lists the SHA-256 of each file and of the
upstream inputs used to build them.

## LINCS2020 provenance clarification (2026-10-25)

The input is the Broad Institute's [public Expanded CMap LINCS Resource 2020](https://clue.io/data/project/5fb7d248c1f885001161d6c8),
specifically [compoundinfo_beta.txt](https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt).
Its [release README](https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/README.txt)
identifies the December 2020 beta release. A fresh public download contained
4,631,014 bytes and matched the frozen input SHA-256
`a71fca6de41dcc46a5063858be7e04155f3c09832c8b3fb35814f03db8d9fdff`.

The [source-specific CLUE redistribution policy](https://clue.io/connectopedia/data_redistribution)
allows redistribution of NIH/publicly funded L1000 data and requests source
and processing attribution. This statement concerns that public release, not
private CLUE resources, access keys, or unrelated licensed content.

L1000PAI normalizes compound identity, gene symbols and action direction into
derived evidence records. It does not bundle this raw metadata or its source
descriptive text. No exact equivalence to a different GEO release is asserted.
The frozen source registry retains its earlier "terms not confirmed" audit
entry for traceability; this dated clarification records the subsequently
verified source and policy. The MIT code license does not override the
separate data terms above.
