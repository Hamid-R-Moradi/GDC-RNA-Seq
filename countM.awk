BEGIN { FS=OFS="\t" }
{ datasets[$1]; fnames[FILENAME]; vals[$1,FILENAME] = $4 }
END {
    printf "%s", "dataset"
    for (fname in fnames) {
        printf "%s%s", OFS, fname
    }
    print ""
    for (dataset in datasets) {
        printf "%s", dataset
        for (fname in fnames) {
            printf "%s%s", OFS, vals[dataset,fname]
        }
        print ""
    }
}
