rule count_individuals_at_step0_from_raw_vcf:
    input:
        vcf = get_vcf_file
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step0.raw.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file"
    threads: 1
    params: 
        step_name = "At step0: raw complete VCF"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step1_selecting_individuals:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step1.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after selecting individuals"
    threads: 1
    params: 
        step_name = "After step1: individual selection"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step2_filtering_on_chromosomes:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.chr.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step2.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after filtering on chromosomes"
    threads: 1
    params: 
        step_name = "After step2: filtering on chromosomes"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step3_biallelic_snps:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step3.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after selecting biallelic SNPs"
    threads: 1
    params: 
        step_name = "After step3: selecting biallelic SNPs"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step4_snp_filters:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step4.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after SNP quality filters"
    threads: 1
    params: 
        step_name = "After step4: SNP quality QC filters"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step5_filter_on_fraction_missing_per_snp:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step5.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after SNP quality filters"
    threads: 1
    params: 
        step_name = "After step5: filtering on fraction missing per SNP"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step6_filter_on_maf:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step6.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after filtering on MAF"
    threads: 1
    params: 
        step_name = "After step6: filtering on MAF"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step7_filter_on_fraction_missing_per_genotype:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step7.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after filtering on fraction missing per genotype"
    threads: 1
    params: 
        step_name = "After step7: filtering on fraction missing per genotype"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

rule count_individuals_after_step8_filter_on_heterozygosity_excess:
    input:
        vcf = WORKING_DIR + "filtered/{sample}.selected.biallelic.qc1.qc2.maf.miss.het.vcf.gz"
    output:
        n_individuals = WORKING_DIR + "counts/{sample}.step8.ind.csv"
    message:
        "Counting number of individuals in {wildcards.sample} VCF file after filtering on heterozygosity excess"
    threads: 1
    params: 
        step_name = "After step8: filtering on heterozygosity excess"
    run:
        count_df = count_individuals(
            vcf_file=input.vcf, 
            n_threads=1, 
            step_name=params.step_name)
        count_df.to_csv(output.n_individuals, index=False)

if config["impute_genotypes"] == "yes":
    rule count_individuals_after_step9_imputation:
        input:
            vcf = RES_DIR + "filtered/{sample}_filtered_imputed.vcf.gz"
        output:
            n_individuals = WORKING_DIR + "counts/{sample}.step9.ind.csv"
        message:
            "Counting number of individuals in {wildcards.sample} VCF file after imputation"
        threads: 1
        params: 
            step_name = "After step9: SNP imputation"
        run:
            count_df = count_individuals(
                vcf_file=input.vcf, 
                n_threads=1, 
                step_name=params.step_name)
            count_df.to_csv(output.n_individuals, index=False)  