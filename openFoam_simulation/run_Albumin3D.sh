#!/bin/bash

set -e

cd "${0%/*}" || exit

echo "=========================================================="
echo " Download and Automatic Configuration"
echo "=========================================================="

# ---------------------------------------------------------
# 1. CONFIGURATION
# ---------------------------------------------------------
#
ZENODO_RECORD_ID="21393950"

ARCHIVE_NAME="Albumin_3D_Case.tar.gz"

# ---------------------------------------------------------
# DOWNLOAD FROM ZENODO
# ---------------------------------------------------------
echo "Download of ${ARCHIVE_NAME}..."
wget -q --show-progress "https://zenodo.org/records/${ZENODO_RECORD_ID}/files/${ARCHIVE_NAME}?download=1" -O "${ARCHIVE_NAME}"

# ---------------------------------------------------------
# UNZIP
# ---------------------------------------------------------
echo "Unzip..."
mkdir -p tmp_zenodo

tar -unzip -q "${ARCHIVE_NAME}" -C tmp_zenodo

# ---------------------------------------------------------
# Files Moving
# ---------------------------------------------------------
echo "Files Moving..."

# Create Folders
mkdir -p simulation_case/constant/triSurface
mkdir -p simulation_case/constant/polyMesh

echo " -> Move file STL in constant/triSurface/..."
find tmp_zenodo -name "*.stl" -exec mv {} simulation_case/constant/triSurface/ \;

echo " -> Mesh Configuration in constant/polyMesh/..."
if find tmp_zenodo -type d -name "polyMesh" | grep -q .; then
    # Move Polymesh Folder
    find tmp_zenodo -type d -name "polyMesh" -exec sh -c 'mv "$1"/* simulation_case/constant/polyMesh/' _ {} \;
fi

echo " -> Initial Conditions Configuration (moving entire 0 folder)..."
if find tmp_zenodo -type d -name "0" | grep -q .; then
    # Move all files from the extracted 0 folder into simulation_case/0/
    find tmp_zenodo -type d -name "0" -exec sh -c 'mv "$1"/* simulation_case/0/' _ {} \;
fi

# ---------------------------------------------------------
# CLEANUP
# ---------------------------------------------------------
echo "Cleanup temporary files..."
rm -rf tmp_zenodo
rm "${ARCHIVE_NAME}"

echo "=========================================================="
echo " Pipeline Complete!"
echo " - STL Files in constant/triSurface/"
echo " - File cell_centers.csv in root Folder"
echo " - Case Ready!"
echo "=========================================================="
