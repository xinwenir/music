#!/bin/bash

# Set the Fortran compiler. It can be modified according to the actual situation, such as changing to 'ifort', etc.
FC=gfortran

# Set the directory where the source files are located and the output directory for the executable file. These can be adjusted according to actual requirements.
SOURCE_DIR=$(dirname $(readlink -f "$0"))/src
echo ${SOURCE_DIR}
OUTPUT_DIR=$(dirname $(readlink -f "$0"))/out

# Define the list of source files, ensuring the file names are accurate and the order conforms to the dependency relationship.
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

# Compilation options. More options can be added as needed, such as optimization levels, etc.
COMPILE_FLAGS="-Wall -Wextra -o ${OUTPUT_DIR}/test-music"

# Check if all source files exist.
for file in "${SOURCE_FILES[@]}"; do
    if [! -f "${SOURCE_DIR}/${file}" ]; then
        echo "Error: The source file ${SOURCE_DIR}/${file} does not exist. Please check the file path and file name."
        exit 1
    fi
done

cd $SOURCE_DIR

# Execute the compilation command.
${FC} ${COMPILE_FLAGS} "${SOURCE_FILES[@]}"

if [ $? -eq 0 ]; then
    echo -e "\033[32mCompilation successful. The executable file ${OUTPUT_DIR}/test-music has been generated.\033[0m"
else
    echo -e "\033[31mCompilation failed. Please check the error message. Possible issues could be syntax errors in the source files or dependency problems, etc.\033[0m"
fi
echo "Running!!!"
cd $OUTPUT_DIR
./test-music