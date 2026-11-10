FROM quay.io/nf-core/ubuntu:20.04

# Install curl (needed for downloading Tracy)
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Download Tracy static binary
RUN curl -L -o /usr/local/bin/tracy https://github.com/gear-genomics/tracy/releases/download/v0.9.1/tracy-v0.9.1-linux-amd64 && \
    chmod +x /usr/local/bin/tracy

# Verify Tracy installation
RUN tracy --version
