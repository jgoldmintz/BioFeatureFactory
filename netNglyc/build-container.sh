#!/bin/bash
#
# Build NetNGlyc Docker Container with ARM64 Mac support
#

set -e

echo "Building NetNGlyc Docker container..."

# Detect if running on ARM64 Mac
if [[ $(uname -m) == "arm64" ]] || [[ $(uname -m) == "aarch64" ]]; then
    echo "Detected ARM64 architecture (Apple Silicon)"
    echo "Container will run in emulation mode"

    # Enable Docker BuildKit for better cross-platform builds
    export DOCKER_BUILDKIT=1
fi

# Clean up any previous extraction
rm -rf netNglyc-1.0 2>/dev/null || true

# Extract fresh NetNGlyc from original tar
echo "Extracting NetNGlyc from original tar..."
tar -xzf "netNglyc.tar.gz"

# Set up NetNGlyc environment
echo "Setting up NetNGlyc environment..."
cd netNglyc-1.0
rm -rf tmp/ 2>/dev/null || true
# Create symlink for Darwin (some scripts check it even on Linux)
ln -s how/how98_Linux how98_Darwin || true
chmod +x how/how98_Linux netNglyc Template/test Template/test_how
cd ..

# Build Docker image with explicit platform
echo "Building Docker image for linux/386 platform..."
docker build --platform linux/386 -t netnglyc:latest .

# Clean up copied files
echo "Cleaning up..."
rm -rf netNglyc-1.0

echo "✅ NetNGlyc Docker container built successfully!"
echo ""
echo "Note: On Apple Silicon Macs, the container will run in emulation mode."
echo "This may be slower but will work correctly."
echo ""
echo "Run with: docker run --rm --platform linux/386 -v \$(pwd):/data netnglyc:latest my_sequences.fasta"