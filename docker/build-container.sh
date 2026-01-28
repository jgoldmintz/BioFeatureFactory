#!/bin/bash
# BioFeatureFactory
# Copyright (C) 2023-2025  Jacob Goldmintz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/bin/bash
#
# Build BioFeatureFactory Software Suite Docker Container with ARM64 Mac support
# Unified container for NetNGlyc, NetPhos, NetMHC, NetSurfP-3.0
#

set -e

echo "Building BioFeatureFactory Software Suite Docker container..."

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

# Set up NetMHC 4.0 environment
echo "Setting up NetMHC 4.0 environment..."
rm -rf netMHC-4.0 2>/dev/null || true
# Copy NetMHC from software directory (2 levels up)
cp -r ../../software/netMHC-4.0 ./netMHC-4.0
if [ -d "netMHC-4.0" ]; then
    echo " NetMHC 4.0 copied successfully"
    chmod +x netMHC-4.0/netMHC
    chmod +x netMHC-4.0/Linux_x86_64/bin/* 2>/dev/null || true
    chmod +x netMHC-4.0/Linux_i686/bin/* 2>/dev/null || true
else
    echo "[!]  Warning: netMHC-4.0 not found in ../../software/ - NetMHC functionality will not work"
fi

# Build Docker image with explicit platform
echo "Building Docker image for linux/386 platform..."
docker build --platform linux/386 -t biofeaturefactory:latest .

# Clean up copied files
echo "Cleaning up..."
rm -rf netNglyc-1.0 netMHC-4.0

echo "BioFeatureFactory Software Suite Docker container built successfully!"
echo ""
echo "Image tagged as: biofeaturefactory:latest"
echo ""
echo "Note: On Apple Silicon Macs, the container will run in emulation mode."
echo "This may be slower but will work correctly."
echo ""
echo "Example usage: docker run --rm --platform linux/386 -v \$(pwd):/data biofeaturefactory:latest my_sequences.fasta"