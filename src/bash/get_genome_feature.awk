function clean_chr(chr, m) {
    if (match(chr, /^chr_([0-9]{1,2})$/, m))
        return "chr" m[1]
    else 
        return "chr"
    }

BEGIN{
    FS = OFS = "\t"
    if (feature_name == "") {
    print "Error: must provide -v utr_type=<pattern>" > "/dev/stderr"
    exit 1
  }
    
}
$0 !~ /^#/ && $1 ~ /^chr_[0-9]{1,2}$/ && tolower($3) ~ tolower(feature_name) {
    $1 = clean_chr($1)
    print $1, $4, $5, $3
    }
