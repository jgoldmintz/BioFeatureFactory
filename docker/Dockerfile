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

# NetNGlyc Docker Container
# Base: 32-bit Linux (i386) for native binary compatibility
FROM --platform=linux/386 debian:bullseye

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install required dependencies including Python for SignalP stub
RUN apt-get update && apt-get install -y \
    perl \
    cpanminus \
    bash \
    gawk \
    grep \
    sed \
    gnuplot \
    netpbm \
    tcsh \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for SignalP stub (using --break-system-packages for Debian 11)
RUN pip3 install --no-cache-dir pandas

# Install required Perl modules
RUN cpanm --notest Array::Base || echo "Array::Base installation failed, continuing"

# Create working directory
WORKDIR /opt/netnglyc

# Copy NetNGlyc distribution (with Darwin symlink already created by build.sh)
COPY netNglyc-1.0/ ./

# Copy SignalP 6.0 stub
COPY signalp_stub /opt/netnglyc/signalp_stub

# Ensure executables have correct permissions
RUN chmod +x netNglyc how/how98_Linux how98_Darwin signalp_stub \
    && chmod +x Template/test Template/test_how

# Environment variables for NetNGlyc
ENV GLYCHOME=/opt/netnglyc
ENV PATH=/opt/netnglyc:$PATH

# Fix netNglyc script paths and configure SignalP with fallback
RUN sed -i 's|^setenv[ \t]*GLYCHOME.*|setenv GLYCHOME /opt/netnglyc|' /opt/netnglyc/netNglyc && \
    # Configure SignalP with intelligent fallback \
    # First, comment out any existing SIGNALP settings \
    sed -i 's|^setenv[ \t]*SIGNALP|#setenv SIGNALP|' /opt/netnglyc/netNglyc && \
    # Add new SignalP configuration with fallback logic after GLYCHOME setting \
    sed -i '/^setenv GLYCHOME/a\\n# SignalP configuration with fallback\n# Try to find signalp6 first, fallback to stub if not found\nif ( -x /usr/local/bin/signalp6 ) then\n   setenv SIGNALP /usr/local/bin/signalp6\nelse if ( -x /usr/bin/signalp6 ) then\n   setenv SIGNALP /usr/bin/signalp6\nelse if ( -x /opt/netnglyc/signalp_stub ) then\n   setenv SIGNALP /opt/netnglyc/signalp_stub\nelse\n   # SignalP not available - NetNGlyc will run without it\n   unsetenv SIGNALP\nendif' /opt/netnglyc/netNglyc

# Create NetNGlyc's required tmp directory with proper permissions
RUN mkdir -p /opt/netnglyc/tmp && chmod 777 /opt/netnglyc/tmp

# Also ensure system temp directory exists with proper permissions
RUN mkdir -p /tmp/netnglyc && chmod 777 /tmp/netnglyc \
    && mkdir -p /data/input /data/output

# Verify SignalP configuration
RUN echo "Checking SignalP configuration in netNglyc script:" && \
    grep -A10 "SignalP configuration" /opt/netnglyc/netNglyc || true && \
    echo "Checking if signalp_stub exists and is executable:" && \
    ls -la /opt/netnglyc/signalp_stub && \
    echo "Testing signalp_stub execution:" && \
    /opt/netnglyc/signalp_stub -h 2>&1 || echo "Stub test completed"

VOLUME ["/data"]

# Entrypoint: run NetNGlyc directly
ENTRYPOINT ["/opt/netnglyc/netNglyc"]
CMD ["-h"]