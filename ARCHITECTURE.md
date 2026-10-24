# L1000PAI repository boundary

This repository is the core R package. It contains package code,
documentation, unit tests, vignettes and small annotation assets only.

The project is split across three sibling workspaces:

| Workspace | Responsibility |
|---|---|
| `L1000PAI` | Core R package and package-level tests |
| `L1000PAIdata` | Downloads and verifies the signature libraries; reference manifest |
| `L1000PAI-analysis` | Case registry, workflows, results manifests and manuscript source |

## Local workspaces

A local checkout may contain large untracked folders, such as historical
analysis outputs. They are excluded from Git and package builds.

## Clean package checks

Windows can encounter path-length failures before `.Rbuildignore` is applied to
large local folders. Use the allowlist staging script instead:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File .\tools\check_paired_packages.ps1 -Quick
```

The script stages only formal package components, builds and installs the
current companion data package into an isolated library, verifies its version,
then builds and checks the core package. It never installs
over the user's normal R library.

Omit `-Quick` for a release-oriented build including vignettes and the manual.
The script leaves its audit directory in the system temporary directory so the
logs remain inspectable.

## Release contract

- Core releases must identify the minimum compatible `L1000PAIdata` version.
- A release candidate must pass unit CI and a paired clean check.
- Package source must not include caches, local checks, logs, analysis outputs
  or unpublished upstream data.
- Signature libraries are never bundled. `L1000PAIdata` downloads them and
  accepts them only when their content matches the reference release.

## Paired-inference boundary

The core package performs directional conjunction testing and can score a
caller-supplied synchronized Up/Down null with `score_joint_activity()`. Null
generation and benchmark orchestration belong to `L1000PAI-analysis` because
they depend on study design, blocking, covariates and source-data provenance.

The core package must not estimate null covariance from observed perturbagen
rows or silently treat Up and Down as independent. When synchronized null draws
are unavailable, strict activity calls use directional conjunction evidence and
the legacy NES difference remains available for ranking.

## Plotting boundary

The plotting functions (`R/plot_*.R`, `R/pathway_enrichment.R`) draw what the
package computes and nothing else: every plot attaches the values it draws as
`attr(p, "plot_data")`, so a figure can be checked against the result table.
Gene-set collections for `pai_enrich()` are kept in the user cache with their
source, version and SHA-256; they are never bundled. Figure layouts made for a
particular paper belong to `L1000PAI-analysis`.

## Report boundary

`pai_deg_report()` (`R/deg_report.R`) composes the package's own scoring,
pathway, annotation and plotting functions into a PDF report and an Excel
workbook; it adds no statistics of its own beyond the checks it reports. It
reads two small reference files from `inst/extdata/report`, documented with
their sources and SHA-256 in the README there and rebuilt by the scripts in
`data-raw/report`: the common-hit percentiles of the Chem library over 12
public GEO perturbation contrasts, and the MSigDB Hallmark E2F targets and G2M checkpoint
gene sets (CC BY 4.0), the only gene sets the package ships.
