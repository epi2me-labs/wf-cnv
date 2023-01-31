# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Changed
- Increased maximum allowed total size of global variables in QDNAseq Rscript

## [v0.0.8]
### Changed
- --genome parameter no longer required (determined from reference or BAM)
- Updated CI command to produce improved report
- Report updated to use histogram function from ezcharts

## [v0.0.7]
### Added
- Workflow can now take BAM as input, as well as FASTQ
### Fixed
- sample_sheet format in schema to expect a file
### Changed
- Documentation update

## [v0.0.6]
### Changed
- Updated description in manifest
- Documentation update

## [v0.0.5]
### Changed
- Documentation update

## [v0.0.4]
### Fixed
- `--sample_sheet` is now a path in the schema
### Changed
- `--map_threads` default is now 8 cores
- Docs update

## [v0.0.3]
### Changed
- QDNAseq R script no longer carries out separate female sample analysis

## [v0.0.2]
### Fixed
- Prevented running with conda profile
### Changed
- Docs update
- nextflow config fix
- Increased memory from 8G -> 16G
- CPU limit for alignment

## [v0.0.1]
### Added
- CNV calling with QDNAseq
