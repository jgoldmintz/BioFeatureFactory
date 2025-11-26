#!/bin/bash
# BioFeatureFactory
# Copyright (C) 2023–2025  Jacob Goldmintz
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