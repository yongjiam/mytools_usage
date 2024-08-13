# Set variables for the file names and wrapping length
input_file="updated_combined.nexus"
output_file="wrapped_updated_combined.nexus"
start_line=6
end_line=27
wrap_length=120

# Extract the lines before the alignment section (up to line 5)
sed -n "1,$((start_line-1))p" $input_file > $output_file

# Wrap the sequence lines from line 6 to 27 and append to the output file
sed -n "${start_line},${end_line}p" $input_file | \
awk -v wrap_len=$wrap_length '{
    name = $1;
    seq = substr($0, length(name) + 2);
    while (length(seq) > wrap_len) {
        print name, substr(seq, 1, wrap_len);
        seq = substr(seq, wrap_len + 1);
        name = "";
    }
    print name, seq
}' >> $output_file

# Extract the lines after the alignment section (starting from line 28)
sed -n "$((end_line+1)),\$p" $input_file >> $output_file
