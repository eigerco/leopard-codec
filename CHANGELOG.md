# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0](https://github.com/eigerco/leopard-codec/releases/tag/leopard-codec-v0.1.0) - 2024-04-15

### Added
- [**breaking**] implement reconstruct for leopard ffe8 ([#5](https://github.com/eigerco/leopard-codec/pull/5))
- leopard8 log and exp lookup tables

### Fixed
- [**breaking**] ensure encoding is compatible with Go's klauspost/reedsolomon ([#4](https://github.com/eigerco/leopard-codec/pull/4))

### Other
- *(ci)* add the release-plz workflow ([#7](https://github.com/eigerco/leopard-codec/pull/7))
- remove Leopard struct and replace it with standalone function ([#3](https://github.com/eigerco/leopard-codec/pull/3))
- add the basic CI configuration ([#2](https://github.com/eigerco/leopard-codec/pull/2))
- remove explicit handling of alignment ([#1](https://github.com/eigerco/leopard-codec/pull/1))
- rename to leopard-codec
- working ffe8 encoding
- wip
- move luts to separate file
- all lookup tables working
- add mul8 lookup table
