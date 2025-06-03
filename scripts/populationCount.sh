#!/bin/bash

awk '
NR > 1 { count[$2]++ }

END {
    print "\\begin{table}[ht]"
    print "\\centering"
    print "\\caption{Population Counts}"
    print "\\label{tab:population_counts}"
    print "\\begin{tabular}{lr}"
    print "\\hline"
    print "Population & Count \\\\"
    print "\\hline"
    for (pop in count) {
        print pop " & " count[pop] " \\\\"}
    print "\\hline"
    print "\\end{tabular}"
    print "\\end{table}"
}' $1 > $2
