#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

directory="./CMake_dlm_hst_plot/FILES"
if [ ! -d "$directory" ]; then
    mkdir "$directory"
    echo "Directory created: $directory"
else
    echo "Directory already exists: $directory"
fi

cd $directory

if ! cmake ../../CMake_dlm_hst_plot; then
	echo "Configuration ${red}failed${reset}"
	return 3
fi

echo "Configuration ${green}successful${reset}"
echo "  To proceed type: make"

cd ../../

return 0
