#!/bin/bash

# Tutorial setup script to prepare the environment for running the tutorial
# Usage: ./tutorial-setup.sh <username>

MyUSERNAME="$1"

set -euo pipefail

mkdir -p ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling
mkdir -p ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling/software
mkdir -p ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling/data

# Download Dorado basecalling models
wget -O ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling/data/models.tar.gz https://example.com/path/to/dorado/models.tar.gz

echo "Tutorial setup completed. Directory structure created and models downloaded."