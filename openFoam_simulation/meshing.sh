#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Ensure we are in the directory where the script is located
cd "${0%/*}" || exit

echo "=========================================================="
echo " Download and Automatic Configuration Manager"
echo "=========================================================="

# ---------------------------------------------------------
# READ USER INPUT AND SELECT CASE
# ---------------------------------------------------------
if [ "$1" == "Albumin" ]; then
    ZENODO_RECORD_ID="20120581"
    ARCHIVE_NAME="Albumin_3D_Case.tar.gz"
elif [ "$1" == "Mapping" ]; then
    ZENODO_RECORD_ID="20120581"
    ARCHIVE_NAME="Dura_Mapping3D_Case.tar.gz"
else
    echo "ERROR: You must specify which dataset to download!"
    echo "Usage:"
    echo "  ./meshing.sh Albumin"
    echo "  or"
    echo "  ./meshing.sh Mapping"
    exit 1
fi

echo "-> Selected case: $1"
echo "-> Archive to download: ${ARCHIVE_NAME}"

# ---------------------------------------------------------
# DOWNLOAD FROM ZENODO
# ---------------------------------------------------------
echo "Downloading..."
wget -q --show-progress "https://zenodo.org/records/${ZENODO_RECORD_ID}/files/${ARCHIVE_NAME}?download=1" -O "${ARCHIVE_NAME}"

# ---------------------------------------------------------
# EXTRACT TO TEMPORARY FOLDER
# ---------------------------------------------------------
echo "Extracting archive..."
mkdir -p tmp_zenodo
tar -xzf "${ARCHIVE_NAME}" -C tmp_zenodo

# ---------------------------------------------------------
# MOVE FILES
# ---------------------------------------------------------
echo "Moving downloaded folders to the local setup..."

# Create local destination folders if they do not exist
mkdir -p simulation_case/constant/triSurface
mkdir -p simulation_case/0

# A. Move constant/triSurface
echo " -> Configuring geometries in constant/triSurface/..."
if find tmp_zenodo -type d -name "triSurface" | grep -q .; then
    # Find the extracted triSurface folder and move its contents to the local one
    find tmp_zenodo -type d -name "triSurface" -exec sh -c 'mv "$1"/* simulation_case/constant/triSurface/' _ {} \;
fi

# B. Move the 0 folder
echo " -> Configuring initial conditions in the 0/ folder..."
if find tmp_zenodo -type d -name "0" | grep -q .; then
    # Find the extracted 0 folder and move its contents to the local one
    find tmp_zenodo -type d -name "0" -exec sh -c 'mv "$1"/* simulation_case/0/' _ {} \;
fi

# ---------------------------------------------------------
# 5. CLEANUP
# ---------------------------------------------------------
echo "Cleaning up temporary files..."
rm -rf tmp_zenodo
rm "${ARCHIVE_NAME}"

echo "=========================================================="
echo " Pipeline Complete!"
echo " The local setup has been successfully updated with the"
echo " files for the $1 case."
echo "=========================================================="
