#!/bin/bash


# Usage:
# ./add_ref_sample.sh input.vcf output.vcf REF_SAMPLE

input_vcf="$1"
output_vcf="$2"
sample_name="$3"

awk -v sample_name="$sample_name" '
BEGIN {OFS="\t"}
/^#CHROM/ {
    $0 = $0 sample_name;
    print;
    next;
}
/^#/ {
    print;
    next;
}
{
    print $0, "0/0:.,.:.:.:.,.,.";
}
' "$input_vcf" > "$output_vcf"
