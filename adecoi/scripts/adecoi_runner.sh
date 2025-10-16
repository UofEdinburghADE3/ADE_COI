

snakemake --snakefile /home/s1680070/repositories/ADE_COI/adecoi/scripts/Snakefile \
            --configfile /home/s1680070/ADE_2025/updated_config.yaml \
            --cores 10 \
            --forceall \
            --rerun-incomplete \
            --nolock

snakemake --snakefile /home/s1680070/repositories/ADE_COI/adecoi/scripts/run_analysis.smk \
            --configfile /home/s1680070/ADE_2025/analysis/analysis_config.yaml \
            --cores 10 \
            --forceall \
            --rerun-incomplete \
            --nolock

# rm /home/s1680070/repositories/ADE_COI/adecoi/ADE_2022/barcode*/consensus_sequences/*.txt

# python adecoi/scripts/index_report.py #need a barcode name map file
