#####################
## SNP counts metrics
#####################

rule count_snps_at_step0_raw_snps:
    input:
        vcf = get_vcf_file
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step0.snp.csv"
    message:
        "Step 0: Counting initial number SNPs in {wildcards.sample} VCF file (all types, SNPs and indels)"
    threads: 1
    params: 
        step_name = "At step 0: N SNPs from raw VCF file"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step1_selecting_individuals:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step1.snp.csv"
    message:
        "Step 1: Counting number of SNPs in {wildcards.sample} VCF file after selecting individuals"
    threads: 1
    params: 
        step_name = "After step 1: N SNPs after individual selection"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step2_filtering_on_chromosomes:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step2.snp.csv"
    message:
        "Step 2: Counting number of SNPs in {wildcards.sample} VCF file after filtering on chromosomes"
    threads: 1
    params: 
        step_name = "After step 2: N SNPs after filtering on chromosomes"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step3_biallelic_snps: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step3.snp.csv"
    message:
        "Step 3: Counting number of biallelic SNPs in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step 3: N biallelic SNPs"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=threads, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step4_snp_filters:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step4.snp.csv"
    message:
        "Step 4: Counting number of SNPs after first filters in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step 4: N SNPs after qc filters"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step5_snp_call_rate_filter: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step5.snp.csv"
    message:
        "Step 5: Counting number of SNPs after filtering on fraction missing per SNP in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step 5: N SNPs after filtering on fraction missing per SNP"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step6_maf_filter_before_imputation: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step6.snp.csv"
    message:
        "Step 6: Counting number of SNPs after MAF filter (before imputation) in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step6: N SNPs after filtering on MAF before imputation"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step7_filter_frac_missing_per_genotype: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step7.snp.csv"
    message:
        "Step 7: Counting number of SNPs after filtering on fraction missing per genotype in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step7: N SNPs after filtering on fraction missing per genotype"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step8_filter_on_heterozygosity_excess: 
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step8.snp.csv"
    message:
        "Step 8: Counting number of SNPs after filtering on heterozygosity excess in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step8: N SNPs after filtering on heterozygosity excess"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step9_imputation: 
    input:
        vcf = RES_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step9.snp.csv"
    message:
        "Step 9: Counting number of SNPs after imputation in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step9: N SNPs after imputation"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)

rule count_snps_after_step10_maf_filter_after_imputation: 
    input:
        vcf = RES_DIR + "filtered/{sample}_filtered_imputed_maf.vcf.gz"
    output:
        n_snps = WORKING_DIR + "counts/{sample}.step10.snp.csv"
    message:
        "Step 10: Counting number of SNPs after MAF filter (after imputation) in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "After step10: N SNPs after filtering on MAF after imputation"
    run:
        count_df = count_variants_by_type(
            vcf_file_path=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_snps, index=False)
