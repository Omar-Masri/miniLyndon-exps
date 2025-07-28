BEGIN {
    LASTR=""
}

(NR % 4 == 2) {
    gsub("-", "", $7)
    R=$7
}

(NR % 4 == 3) {
    if (R != LASTR) {
        print ">" $2
        print R
        LASTR=R
    }
}
