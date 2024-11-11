#!/bin/bash

# Set the Fortran compiler. Here we assume using gfortran, which can be modified according to the actual installation situation.
FC=gfortran

# Set the directory where the source files are located. Here we assume all source files are in the directory named'music', which can be adjusted as needed.
SOURCE_DIR="src"

# Define the list of source files, listing each source file in the order of dependency relationships.
SOURCE_FILES=(
    "test-music.f"
    "music-crosssections.f"
    "ranlux.f"
    "music.f"
    "corset.f"
    "corgen.f"
    "rnormal.f"
    "ranmar.f"
)

# Compilation options. Here we add some common options, such as displaying warning messages, etc. These can be further expanded according to requirements.
COMPILE_FLAGS="-Wall -Wextra -o test-music"

# Check if all source files exist in the specified directory.
for file in "${SOURCE_FILES[@]}"; do
    if [ ! -f "${SOURCE_DIR}/${file}" ]; then
        echo "Error: The source file ${SOURCE_DIR}/${file} does not exist. Please check the file path and file name."
        exit 1
    fi
done

# Execute the compilation command to compile all source files according to the set compilation options.
${FC} ${COMPILE_FLAGS} "${SOURCE_DIR}/${SOURCE_FILES[@]}"

# Determine whether the compilation is successful based on the return value of the compilation command.
if [ $? -eq 0 ]; then
    echo "Compilation successful. The executable file test-music has been generated."
else
    echo "Compilation failed. Please check the error message. Possible issues could be syntax errors in the source files or dependency problems, etc."
fi
