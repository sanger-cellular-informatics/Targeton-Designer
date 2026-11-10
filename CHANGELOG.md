# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project follows Semantic Versioning where applicable.

---
## [v1.1.0] - November 10th, 2026

### Added
- Option to retrieve sequences from a local reference genome, the path to which is specified with `reference_genome`.
FASTA headers with or without a `chr` prefix (e.g. `1` or `chr1`) are supported

### Changed
- Updated to use v0.0.4 of `cicd-templates` and runner tag updated
- Updated the Python Docker base image OS from Debian 11 (Bullseye) to Debian 12 (Bookworm)



## [v1.0.0] - April 14, 2026

### Added
- Output primer pairs in ipcress-compatible format when "ipcress_parameters" are set up in designer configuration.
- Use kmer-list masking to prevent designing primers in repetitive or non-specific regions.
- Rank primers based on kmer-masking stringency (using the PRIMER_MASK_FAILURE_RATE Primer3 parameter).
- Support for ranking primers by product size.
- Implemented additional filtering capabilities for:
  - duplicates
  - HAP1_variant

  These filters help refine primer selection and improve the quality of generated results.

### Changed
- Improved primer selection workflow with enhanced ranking and filtering logic.

### Removed/Deprecated
The following commands are no longer in use:
- scoring
- slicer
- generate_targeton_csv
- collate_primer_data
- post_primers

---

## [v0.1.0] - June 28, 2022

### Added
- Initial release for UAT.
- Delivered standalone Slicer and Primer3 tools as standalone modules.
- Introduced base functionality for running tools in a UAT environment.

### Notes
- Sprint 7 increment.
- Initial release enabling UAT runs through the delivered modules.
- Evidence collection file: `v0.1.0-evidences-150.json` (collected 28 Jun 2022).

---
