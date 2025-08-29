# Docker-based NetNGlyc for Disease-Associated Synonymous Mutations

This Docker setup provides accurate NetNGlyc predictions with SignalP 6.0 integration for analyzing N-glycosylation sites in disease-associated synonymous mutations.

## 🍎 Apple Silicon Solution

**Primary Use Case**: This Docker solution is specifically designed to solve NetNGlyc compatibility issues on **Apple Silicon (M1/M2/M3) Macs**. 

- **Problem**: NetNGlyc's 32-bit Linux binaries cannot run natively on ARM64 macOS
- **Solution**: Docker provides a 32-bit Linux environment where NetNGlyc works perfectly
- **Cross-platform**: Also works on Intel Macs and Linux systems for consistency

## 📋 Required Files

### Essential Files (Required)
- **`Dockerfile`** - Container definition with 32-bit Linux environment
- **`build-container.sh`** - Automated build script with ARM64 Mac support
- **`signalp_stub`** - SignalP 6.0 integration script
- **`full-docker-netnglyc-pipeline.py`** - Complete production pipeline
- **`netNglyc.tar.gz`** - Original NetNGlyc 1.0 distribution (download from DTU)

### Optional Files (Recommended)
- **`docker_netnglyc_parallel.py`** - Parallel processing module
- **`pipeline_docker_integration.py`** - Integration bridge for existing workflows

## 🚀 Quick Start

### 1. Obtain Required Software

#### NetNGlyc 1.0
Download `netNglyc.tar.gz` from:
- **Official source**: https://services.healthtech.dtu.dk/software.php
- **License**: Academic use only - requires institutional email registration
- **Place in**: `/Volumes/990pro/disease_associated_synonymous_mutations/docker/`

#### SignalP 6.0 (Required for signalp_stub)
-  Download and install SignalP 6.0 from:
- **Official source**: https://services.healthtech.dtu.dk/software.php
- **License**: Academic use only - requires institutional email registration

- **Install locatialson**: Standard system path (e.g., `/usr/local/bin/signalp6` or conda environment)
- **Note**: The `signalp_stub` script calls SignalP 6.0 internally for signal peptide predictions

### 2. Build Docker Container
```bash
cd /Volumes/990pro/disease_associated_synonymous_mutations/docker
./build-container.sh
```

### 3. Run Complete Pipeline
```bash
# Process both wildtype and mutant sequences
python3 full-docker-netnglyc-pipeline.py \
    --fasta_wt ../FASTA_files/wt/aaseq \
    --fasta_mut ../FASTA_files/mut/aaseq \
    --output_wt ../Data/netnglyc-docker-out-wt \
    --output_mut ../Data/netnglyc-docker-out-mut \
    --mapping_dir ../mutations/combined/aa \
    --workers 4
```

## 🛠️ Configuration Options

### Adjust Worker Count
```bash
# Faster processing (more RAM usage)
python3 docker_netnglyc_parallel.py INPUT OUTPUT --workers 8

# Slower but less resource intensive
python3 docker_netnglyc_parallel.py INPUT OUTPUT --workers 2
```

### Process Specific Files
```bash
# Only process specific FASTA pattern
python3 docker_netnglyc_parallel.py INPUT OUTPUT --pattern "ABCB1*.fasta"
```

## 🧪 Testing & Validation

### Quick Test
```bash
python3 test_single_sequence.py
```

### Validate Against Known Results
```bash
# Compare with your original 164 predictions
# The Docker version should find the same glycosylation sites
```

## 📁 File Structure

```
docker/
├── Dockerfile                           # Container definition (32-bit Linux)
├── build-container.sh                   # Automated build script
├── signalp_stub                         # SignalP 6.0 integration
├── full-docker-netnglyc-pipeline.py     # Complete production pipeline
├── docker_netnglyc_parallel.py          # Parallel processing module
├── pipeline_docker_integration.py       # Integration bridge
├── netNglyc.tar.gz                      # NetNGlyc 1.0 distribution (download)
└── README.md                            # This documentation
```

## 🔧 Integration Options

### Option 1: Standalone Docker Pipeline
```bash
# Complete self-contained processing
python3 full-docker-netnglyc-pipeline.py --help
```

### Option 2: Integration with Existing Pipeline
```python
# Add to existing run-netnglyc-pipeline.py
import sys
sys.path.append('docker')
from pipeline_docker_integration import run_netnglyc_sequential_docker
```

## ⚠️ Prerequisites

1. **Docker Desktop**: Install from https://www.docker.com/products/docker-desktop
2. **NetNGlyc License**: Academic registration required at DTU website  
3. **SignalP 6.0**: Must be installed and accessible in system PATH
4. **4-8GB RAM**: For parallel processing (adjustable)
5. **Storage**: ~2GB for container + temporary processing files

## 🐛 Troubleshooting

### Docker not found
```bash
# Install Docker Desktop and ensure it's running
docker --version
```

### Container build fails
```bash
# Check NetNGlyc source exists
ls -la ../software/netNglyc-1.0/
```

### Processing fails
```bash
# Test container setup
python3 docker_netnglyc_parallel.py --test
```

### Slow performance
```bash
# Reduce workers if running out of memory
python3 docker_netnglyc_parallel.py INPUT OUTPUT --workers 2
```

## 📈 Expected Results

With this Docker setup, you should get:
- ✅ **Identical predictions** to your original 164 results
- ✅ **Working neural network** (no more "Command not found" errors)
- ✅ **SignalP 6.0 integration** with accurate signal peptide predictions
- ✅ **Parallel processing** for faster completion

The Docker approach ensures 100% compatibility with the original NetNGlyc while providing modern SignalP 6.0 integration.