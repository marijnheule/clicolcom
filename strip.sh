cat - | tr " " "\n" | awk '{if ($1 > 0) print $0}' | grep -v "v"
#cat $1 | tr " " "\n" | awk '{if ($1 > 0) print $0}' | grep -v "v"
