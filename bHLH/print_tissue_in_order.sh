head -n1 *.csv |grep "^\""|sed 's/,/|/1'|cut -d'|' -f2|awk -F, '{for(i=1; i<=NF; i++) {split($i, a, "."); if (a[2]) printf "%s ", a[2]} print ""}'
