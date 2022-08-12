#!/bin/bash
# get user input with prompt to the path where the zipped fmri data lies 
echo "Please enter the full path to the zipped meg data file:"
read -e -p "Path: " data_path
echo "You entered: "$data_path""

# get user input with prompt to the path where the cloned github directory lies
echo "Please enter the path to the cloned github directory:"
read -e -p "Path: " github_path
echo "You entered: $github_path"
# create a folder inside the cloned github directory and unzip the files to that folder
echo "Creating a folder inside the cloned github directory and unzipping files to that folder"
mkdir -p $github_path/data/meg/
# unzip the files in the given path 
echo "Unzipping the files in the given path"
unzip -q "$data_path" -d $github_path/data/meg/
echo "Unzipping the single subject MEG data"
cd $github_path/data/meg/preproc/
tar -xvf $github_path/data/meg/preproc/meg_preproc.tar.gz