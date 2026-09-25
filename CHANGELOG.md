# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-30
### Added
- Add the first workable version.

## [1.0.3] - 2026-01-31
### Added
- Added annotate sage vcf with purple

## [1.0.4] - 2026-02-10
- update sage version to 4.2, wisp version to 1.2, so sage vcf will
includes AED will field which will be a filter in wisp
- sage_append use purple vcf as input instead of sage_primary
- removed annotate plasma with purple task
- removed merge bqr task
- wisp use skip_bqr option, and will provision out somatic_variant_tsv
separately

## [1.0.5] - 2026-02-12
### Added
- Added PAVE task

### Changed
- SAGE append not run in parallel
- WISP uses bqr_dor from sage_append

## [3.0.0] - 2026-09-16
### Added
- A brand new version based on oncoanalyser/3.0.0

## [3.1.0] - 2026-09-24
### Added
- Label each row of the summary with the sample it came from.

### Fixed
- Pass the sequencing platform to WISP, which otherwise applies its Illumina
error model to every sample whatever the run.

## [3.2.0] - 2026-09-25
### Changed
- Measure every sample in one SAGE append and one WISP call rather than a pair per
sample, which the tool supports once the sample list is separated by semicolons.
- Refuse a run that asks for copy-number evidence alongside further samples.
