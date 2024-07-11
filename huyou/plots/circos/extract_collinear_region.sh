#!/bin/bash

# Define the input file
input_file="hap1.collinearity"

# Initialize variables
group_start="## Alignment"
first_line=""
last_line=""

# Read the file line by line
while IFS= read -r line; do
    if [[ $line == $group_start* ]]; then
        # Print the first and last lines of the previous group if they exist
        if [[ -n $first_line ]]; then
            echo "$first_line"
            echo "$last_line"
        fi
        # Reset variables for the new group
        first_line=""
        last_line=""
    else
        # If first_line is empty, this is the first data line
        if [[ -z $first_line ]]; then
            first_line="$line"
        fi
        # Always update the last_line
        last_line="$line"
    fi
done < "$input_file"

# Print the first and last lines of the last group
if [[ -n $first_line ]]; then
    echo "$first_line"
    echo "$last_line"
fi
