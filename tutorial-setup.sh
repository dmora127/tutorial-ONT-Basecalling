#!/bin/bash

# Tutorial setup script to prepare the environment for running the tutorial
# Usage: ./tutorial-setup.sh <username>

MyUSERNAME="$1"

set -euo pipefail

mkdir -p /ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling
mkdir -p /ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling/software
mkdir -p /ospool/ap40/data/$MyUSERNAME/tutorial-ONT-Basecalling/data

echo "Tutorial setup completed. Directory structure created."