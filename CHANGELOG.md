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