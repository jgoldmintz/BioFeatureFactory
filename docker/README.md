# Shared Docker Infrastructure

This directory contains the shared Docker infrastructure used by both NetPhos and NetNGlyc prediction pipelines.

## Files

- **Dockerfile** - Container definition with 32-bit Linux environment for running NetNGlyc and NetPhos
- **build-container.sh** - Automated build script with ARM64 Mac support
- **signalp_stub** - SignalP 6.0 integration script for signal peptide predictions
- **netnglyc_stub** - License-free NetNGlyc alternative for fallback scenarios

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

The container provides a 32-bit Linux environment that enables:
- **Modern SignalP 6.0 integration** for accurate signal peptide predictions
- **Apple Silicon compatibility** where native NetNGlyc and NetPhos binaries cannot run
- **Consistent cross-platform execution** across different architectures

## Platform Requirements

### **Docker Required:**
- **Apple Silicon Macs** (M1/M2/M3) - Native binaries incompatible with ARM64
- **Windows** - Linux binary compatibility

### **Docker Optional:**
- **Linux x86_64** - Can run NetNGlyc natively with 32-bit libraries
- **Intel Macs** - Can run Darwin binaries natively (if available)

### **Prerequisites:**
1. **SignalP 6.0** Install on a host system for modern signal peptide prediction
2. **Docker Desktop** (Platform-dependent) - Required for Apple Silicon and Windows
3. **NetNGlyc 1.0** - Download and place netNglyc.tar.gz in this directory

## Container Image

The built container is tagged as netnglyc:latest and contains:
- NetNGlyc 1.0 prediction software
- APE (A Package of Programs for Sequence Analysis) 
- 32-bit Linux environment compatibility
- Integration with host SignalP 6.0 via signalp_stub