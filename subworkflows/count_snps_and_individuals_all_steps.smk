############
# SNP counts
############

if config["impute_genotypes"] == "yes":
    rule merge_snp_counts_from_all_steps: 
        input:
            expand(WORKING_DIR + "counts/{sample}.step{step}.snp.csv",
            sample=SAMPLES, 
            step=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        output:
            RES_DIR + "counts/counts_merged.snp.csv"
        message:
            "Merging all SNP counts from different steps into a summary file"
        params: 
            out_path = RES_DIR + "counts/counts_merged.snp.csv"
        threads: 1
        run:
            # create an empty dataframe with Pandas and add each count dataframe to it
            counts_df = pd.DataFrame()
            for f in input:
                df = pd.read_csv(f, index_col=0).head() 
                counts_df = pd.concat([counts_df, df], axis=0)
            print(counts_df)
            counts_df.to_csv(path_or_buf=params.out_path, index=True)
elif config["impute_genotypes"] == "no":
    rule merge_snp_counts_from_all_steps: 
        input:
            expand(WORKING_DIR + "counts/{sample}.step{step}.snp.csv",
            sample=SAMPLES, 
            step=[0, 1, 2, 3, 4, 5, 6, 7, 8])
        output:
            RES_DIR + "counts/counts_merged.snp.csv"
        message:
            "Merging all SNP counts from different steps into a summary file"
        params: 
            out_path = RES_DIR + "counts/counts_merged.snp.csv"
        threads: 1
        run:
            # create an empty dataframe with Pandas and add each count dataframe to it
            counts_df = pd.DataFrame()
            for f in input:
                df = pd.read_csv(f, index_col=0).head() 
                counts_df = pd.concat([counts_df, df], axis=0)
            counts_df.to_csv(path_or_buf=params.out_path, index=True)

###################
# Individual counts
###################

if config["impute_genotypes"] == "yes":
    rule merge_individual_counts_from_all_steps:
        input:
            expand(WORKING_DIR + "counts/{sample}.step{step}.ind.csv",
                   sample=SAMPLES,
                   step=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        output:
            RES_DIR + "counts/counts_merged.ind.csv"
        message:
            "Merging all individuals counts from different steps into a summary file"
        params:
            out_path = RES_DIR + "counts/counts_merged.ind.csv"
        threads: 1
        run:
            counts_df = pd.DataFrame()
            for f in input:
                df = pd.read_csv(f, index_col=0)
                # add lines to df
                counts_df = pd.concat([counts_df, df], axis=0)
            counts_df.to_csv(path_or_buf=params.out_path, index=True)
elif config["impute_genotypes"] == "no":
    rule merge_individual_counts_from_all_steps:
        input:
            expand(WORKING_DIR + "counts/{sample}.step{step}.ind.csv",
                   sample=SAMPLES,
                   step=[0, 1, 2, 3, 4, 5, 6, 7, 8])
        output:
            RES_DIR + "counts/counts_merged.ind.csv"
        message:
            "Merging all individuals counts from different steps into a summary file"
        params:
            out_path = RES_DIR + "counts/counts_merged.ind.csv"
        threads: 1
        run:
            counts_df = pd.DataFrame()
            for f in input:
                df = pd.read_csv(f, index_col=0)
                # add lines to df
                counts_df = pd.concat([counts_df, df], axis=0)
            counts_df.to_csv(path_or_buf=params.out_path, index=True)
