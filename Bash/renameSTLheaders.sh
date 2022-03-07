#!/bin/bash

find ./ -name \*.stl -type f -print0 | while read -d $'\0' file; do
    echo "Processing $file"
    fileName="$(basename -s .stl $file)"
    sed -i "s/Exported from Blender-2.90.1/$fileName/g" $file
done