for file in *.csv; do echo "$file: $(head -n 1 "$file" | awk -F, '{print NF}') columns"; done
