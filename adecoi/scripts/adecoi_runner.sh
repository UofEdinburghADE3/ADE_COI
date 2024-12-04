

snakemake --snakefile /localdisk/home/s1680070/repositories/ADE_COI/adecoi/scripts/Snakefile \
            --configfile /localdisk/home/s1680070/repositories/ADE_COI/adecoi/data/config.yaml \
            --cores 10 \
            --forceall \
            --rerun-incomplete \
            --nolock

snakemake --snakefile /localdisk/home/s1680070/repositories/ADE_COI/adecoi/scripts/run_analysis.smk \
            --configfile /localdisk/home/s1680070/ADE_2022/combined_analysis/analysis_config.yaml \
            --cores 10 \
            --forceall \
            --rerun-incomplete \
            --nolock

rm /localdisk/home/s1680070/repositories/ADE_COI/adecoi/ADE_2022/barcode*/consensus_sequences/*.txt

python adecoi/scripts/index_report.py #need a barcode name map file