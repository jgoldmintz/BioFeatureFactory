# Shared Docker Infrastructure

This directory contains the shared Docker infrastructure used by both NetPhos and NetNGlyc prediction pipelines.

## Files

- **Dockerfile** - Container definition with 32-bit Linux environment for running NetNGlyc and NetPhos
- **build-container.sh** - Automated build script with ARM64 Mac support
- **signalp_stub** - SignalP 6.0 integration script for signal peptide predictions
- **netNglyc_stub** - NetNglyc 1.0 integration script for N-linked glycolysation
## Required Downloads

### NetNGlyc 1.0 Distribution
Download netNglyc.tar.gz and place it in this directory:
- **Source**: https://services.healthtech.dtu.dk/software.php
- **License**: Academic use only - requires institutional email registration
- **Placement**: Place the netNglyc.tar.gz file in this docker/ directory

## Building the Container

```bash
cd path/to/BioFeatureFactory/docker/
./build-container.sh
```

## Container Usage

This Docker container is used by:
- **NetNGlyc Pipeline** (../netnglyc/full-docker-netnglyc-pipeline.py)
- **NetPhos Pipeline** (../netphos/netphos-pipeline.py)

The container provides a 32-bit Linux environment where both NetNGlyc and NetPhos can run correctly, especially important for Apple Silicon Macs where the original binaries cannot run natively.

## Prerequisites

1. **Docker Desktop**: Install from https://www.docker.com/products/docker-desktop
2. **SignalP 6.0**: Install on host system for signalp_stub integration
3. **NetNGlyc 1.0**: Download and place netNglyc.tar.gz in this directory

## Container Image

The built container is tagged as netnglyc:latest and contains:
- NetNGlyc 1.0 prediction software
- APE (A Package of Programs for Sequence Analysis) 
- 32-bit Linux environment compatibility
- Integration with host SignalP 6.0 via signalp_stub
