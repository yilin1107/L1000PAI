# Frozen catalog byte format

Plain CSV catalogs use UTF-8 without a byte-order mark and LF line endings,
consistent with this repository's `.gitattributes`. Manifest SHA-256 values
refer to the exact shipped bytes, not a platform-normalized interpretation.
Compressed CSV files remain binary and are not reserialized.

For the 1.1.0 release, the four v2 CSV tables were normalized from CRLF to LF.
No field values, row order, taxonomy memberships or scientific results changed.
The v2 manifest and the shared v2 identity-table entry in the v3 manifest now
record the LF hashes. The original CRLF hashes are retained below for provenance:

| Artifact | Original CRLF SHA-256 |
|---|---|
| chem_term_registry_v2 | 26dc805d6ffae6f6aebabda80bfcbecc28b75c78a55fa98689564c0e0d6d06e0 |
| chem_family_registry_v2 | 5fb530c5aabe69a0094c7a246f3765d87185e04b822b558623bdda39f243358c |
| chem_family_membership_v2 | 2abc6eafb91ed6c212805faab059f27183e15b096d92233b6cfb754b44f2c0c4 |
| chem_source_registry_v2 | 89faf12e30d1ff602d8d3506926d95f5b49da91d52a614261f71189341b96a00 |

The manifest's creation dates describe the original scientific tables; this
packaging-only normalization was verified on 2026-10-25. Strict byte-hash
verification remains enabled, including on the shared identity table.
