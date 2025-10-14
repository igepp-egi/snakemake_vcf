if config["impute_genotypes"] == "yes":
    rule extract_individuals_and_genotypes_from_final_vcf:
        input: 
            vcf = RES_DIR + "filtered/{sample}_filtered_imputed_maf.vcf.gz"
        output: 
            individuals = WORKING_DIR + "genotypes_final_step/{sample}.list_individuals.csv",
            genotypes = WORKING_DIR + "genotypes_final_step/{sample}.genotypes.GT.FORMAT" 
        message:
            "Extracting genotypes and individuals for {wildcards.sample} from {input.vcf} VCF file"
        params: 
            out_prefix = WORKING_DIR + "genotypes_final_step/{sample}.genotypes"
        threads: 1
        shell:
            "bcftools query --list-samples {input.vcf} > {output.individuals};"
            "vcftools --gzvcf {input.vcf} --extract-FORMAT-info GT --out {params.out_prefix} "
elif config["impute_genotypes"] == "no":
    rule extract_individuals_and_genotypes_from_final_vcf:
        input: 
            vcf = RES_DIR + "filtered/{sample}_filtered_not_imputed.vcf.gz"
        output: 
            individuals = WORKING_DIR + "genotypes_final_step/{sample}.list_individuals.csv",
            genotypes = WORKING_DIR + "genotypes_final_step/{sample}.genotypes.GT.FORMAT" 
        message:
            "Extracting genotypes and individuals for {wildcards.sample} from {input.vcf} VCF file"
        params: 
            out_prefix = WORKING_DIR + "genotypes/{sample}.genotypes"
        threads: 1
        shell:
            "bcftools query --list-samples {input.vcf} > {output.individuals};"
            "vcftools --gzvcf {input.vcf} --extract-FORMAT-info GT --out {params.out_prefix} "

rule create_final_genotypes_tsv:
    input:
        individuals = WORKING_DIR + "genotypes_final_step/{sample}.list_individuals.csv",
        genotypes = WORKING_DIR + "genotypes_final_step/{sample}.genotypes.GT.FORMAT"
    output:
        RES_DIR + "genotypes/{sample}.genotypes.tsv"
    message:
        "Create the {wildcards.sample} genotype tsv file from the final filtered VCF"
    params:
        output_file_path = RES_DIR + "genotypes/{sample}.genotypes.tsv"
    threads: 1
    run:
        create_genotype_tsv(
            genotypes_tsv=input.genotypes, 
            individuals_csv=input.individuals, 
            output_tsv=params.output_file_path)