## replace the ID with names in the predicted fasta
awk '/^>/ {for (i=1; i<=NF; i++) if ($i ~ /Name=/) print ">"$i; next} {print}' predicted_proteins.fasta > simple_predicted_proteins.fasta
