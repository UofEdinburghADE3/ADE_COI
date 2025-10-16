import os
from Bio import SeqIO
from adecoi.utils import histograms
from adecoi.utils import kraken_parsing
from adecoi.utils import report_runner as report
import gzip
import collections
import yaml


rule all:
    input:
        expand(os.path.join(config["output_path"],"{barcode}","classified","krona.filtered.html"), barcode=config["barcodes"]),
        expand(os.path.join(config["output_path"],"{barcode}","classified","krona.unfiltered.html"), barcode=config["barcodes"]),
        expand(os.path.join(config["output_path"],"{barcode}","processed_taxa","taxa_present.csv"), barcode=config["barcodes"]),
        expand(os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt"), barcode=config["barcodes"]),
        expand(os.path.join(config["output_path"],"{barcode}","analysis_report.html"), barcode=config["barcodes"]),
        expand(os.path.join(config["repo_path"],"{barcode}","krona.unfiltered.html"), barcode=config["barcodes"]),
        expand(os.path.join(config["repo_path"],"{barcode}","consensus_sequences","cns.prompt.txt"),barcode=config["barcodes"]),
        expand(os.path.join(config["publish_path"],"{barcode}","consensus_sequences","cns.prompt.txt"), barcode=config["barcodes"])


rule kraken_classify:
    input:
        db=config["kraken_db"],
        binned=os.path.join(config["output_path"],"{barcode}","unfiltered","nanopore_reads.fastq"),
        taxonomy="/home/s1680070/seawater_wimp/kraken_data/ADE_working/classified/taxonomy.tab"
    output:
        kraken=os.path.join(config["output_path"],"{barcode}","classified","unfiltered.kraken"),
        classified_out = os.path.join(config["output_path"],"{barcode}","classified","classified_reads.unfiltered.fastq"),
        report = os.path.join(config["output_path"],"{barcode}","classified","kraken_report.unfiltered.kreport2")
    params:
        barcode="{barcode}",
        classified=os.path.join(config["output_path"],"{barcode}","classified")
    threads:
        2
    shell:
        """
        kraken2 --db {input.db} \
            --report {output.report} \
            --classified-out {output.classified_out} \
            --memory-mapping \
            {input.binned} \
            > {output.kraken}
        """

rule make_krona:
    input:
        kraken=os.path.join(config["output_path"],"{barcode}","classified","unfiltered.kraken")
    output:
        krona=os.path.join(config["output_path"],"{barcode}","classified","krona.unfiltered.html")
    log: os.path.join(config["output_path"],"{barcode}","classified","krona.unfiltered.log")
    shell:
        """
        ktImportTaxonomy \
            -q 2 -t 3 {input.kraken} \
            -o {output.krona} \
            &> {log} 
        """

rule filter_by_length:
    input:
        os.path.join(config["output_path"],"{barcode}","unfiltered","nanopore_reads.fastq")
    output:
        reads = os.path.join(config["output_path"],"{barcode}","read_length_filtered","filtered_reads.fastq"),
        histogram1 = os.path.join(config["output_path"],"{barcode}","figures","unfiltered_reads_histogram.svg"),
        histogram2 = os.path.join(config["output_path"],"{barcode}","figures","length_filtered_reads_histogram.svg")
    run:
        fastq_records = []
        lengths = []
        filtered_lengths = []
        with open(output.reads,"w") as fw:
            for record in SeqIO.parse(input[0],"fastq"):
                length = len(record)
                lengths.append(length)
                if length > int(config["min_length"]) and length < int(config["max_length"]):
                    fastq_records.append(record)
                    filtered_lengths.append(length)
            SeqIO.write(fastq_records,fw, "fastq")
        histograms.make_read_length_histogram_svg(lengths,config["min_length"],config["max_length"], output.histogram1)
        histograms.make_read_length_histogram_svg(filtered_lengths,config["min_length"],config["max_length"], output.histogram2)

rule kraken_classify_filtered:
    input:
        db=config["kraken_db"],
        binned=os.path.join(config["output_path"],"{barcode}","read_length_filtered","filtered_reads.fastq"),
        taxonomy="/home/s1680070/seawater_wimp/kraken_data/ADE_working/classified/taxonomy.tab"
    output:
        kraken=os.path.join(config["output_path"],"{barcode}","classified","filtered.kraken"),
        classified_out = os.path.join(config["output_path"],"{barcode}","classified","classified_reads.filtered.fastq"),
        report = os.path.join(config["output_path"],"{barcode}","classified","kraken_report.filtered.kreport2")
    params:
        barcode="{barcode}",
    threads:
        2
    shell:
        """
        kraken2 --db {input.db} \
            --report {output.report} \
            --classified-out {output.classified_out} \
            --memory-mapping \
            {input.binned} \
            > {output.kraken}
        """

rule make_krona_filtered:
    input:
        kraken=os.path.join(config["output_path"],"{barcode}","classified","filtered.kraken")
    output:
        krona=os.path.join(config["output_path"],"{barcode}","classified","krona.filtered.html")
    log: os.path.join(config["output_path"],"{barcode}","classified","krona.filtered.log")
    shell:
        """
        ktImportTaxonomy \
            -q 2 -t 3 {input.kraken} \
            -o {output.krona} \
            &> {log} 
        """

rule get_taxa_reads:
    input:
        report = os.path.join(config["output_path"],"{barcode}","classified","kraken_report.filtered.kreport2"),
        classified_out = os.path.join(config["output_path"],"{barcode}","classified","classified_reads.filtered.fastq")
    params:
        barcode = "{barcode}"
    output:
        report = os.path.join(config["output_path"],"{barcode}","processed_taxa","taxa_present.csv"),
        yaml = os.path.join(config["output_path"],"{barcode}","config.yaml")
    run:
        taxa_dict = kraken_parsing.parse_report_for_taxa(input.report,config["min_reads"],output.report)
        gathered_taxa_read_path = os.path.join(config["output_path"],f"{params.barcode}","taxa_read_files")
        if not os.path.exists(gathered_taxa_read_path):
            os.mkdir(gathered_taxa_read_path)
        reference_path = os.path.join(config["output_path"],f"{params.barcode}","references")
        if not os.path.exists(reference_path):
            os.mkdir(reference_path)
        config["reference_path"] = reference_path
        config["read_path"] = gathered_taxa_read_path
        taxa = kraken_parsing.get_reads(taxa_dict,input.classified_out,gathered_taxa_read_path)
        database = config["kraken_fasta"]
        map_file = os.path.join(config["kraken_db"],"seqid2taxid.map")

        valid_taxa = kraken_parsing.extract_taxid_ref(taxa,map_file,database,reference_path)
        config["taxa"] = valid_taxa
        with open(output.yaml, 'w') as fw:
            yaml.dump(config, fw) 

rule get_consensus_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"cns_runner.smk"),
        taxa = os.path.join(config["output_path"],"{barcode}","processed_taxa","taxa_present.csv"),
        yaml = os.path.join(config["output_path"],"{barcode}","config.yaml")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences")
    threads: 4
    output:
        taxa = os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt")
    run:
        print("Running consensus pipeline.")
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} "
                    "--cores 4 && touch {output.taxa}")

rule generate_report:
    input:
        all_reads = os.path.join(config["output_path"],"{barcode}","unfiltered","nanopore_reads.fastq"),
        classified_reads = os.path.join(config["output_path"],"{barcode}","classified","classified_reads.unfiltered.fastq"),
        filtered_reads = os.path.join(config["output_path"],"{barcode}","read_length_filtered","filtered_reads.fastq"),
        classified_filtered = os.path.join(config["output_path"],"{barcode}","classified","classified_reads.filtered.fastq"),
        taxa = os.path.join(config["output_path"],"{barcode}","processed_taxa","taxa_present.csv"),
        prompt = os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt"),
        yaml = os.path.join(config["output_path"],"{barcode}","config.yaml"),
        histogram1 = os.path.join(config["output_path"],"{barcode}","figures","unfiltered_reads_histogram.svg"),
        histogram2 = os.path.join(config["output_path"],"{barcode}","figures","length_filtered_reads_histogram.svg")
    params:
        barcode = "{barcode}",
        cns_path = os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences")
    output:
        html = os.path.join(config["output_path"],"{barcode}","analysis_report.html")
    run:
        # for seq in os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences"):
        # if file empty, remove, else add to table
        data_for_report = {}
        report.get_read_counts(input.all_reads, input.classified_reads, input.filtered_reads, input.classified_filtered,data_for_report)
        report.load_histograms_svg(input.histogram1, input.histogram2, data_for_report)

        to_remove = report.data_for_table(input.taxa,params.cns_path, data_for_report)

        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)
        try:
            year = config["year"]
        except:
            print("not in config")
            year = "2022"
        report.make_report(output.html,config_loaded,data_for_report,params.barcode,year)
        # for i in to_remove:
        #     shell(f"rm {i}")



rule publish:
    input:
        report = os.path.join(config["output_path"],"{barcode}","analysis_report.html"),
        fkrona = os.path.join(config["output_path"],"{barcode}","classified","krona.filtered.html"),
        ukrona = os.path.join(config["output_path"],"{barcode}","classified","krona.unfiltered.html"),
        cns= os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt")
    params:
        barcode = "{barcode}",
        cns_path = os.path.join(config["output_path"],"{barcode}","processed_taxa","consensus_sequences","consensus_sequences"),
        rcns_path = os.path.join(config["repo_path"],"{barcode}","consensus_sequences"),
        pcns_path = os.path.join(config["publish_path"],"{barcode}","consensus_sequences")
    output:
        preport = os.path.join(config["publish_path"],"{barcode}","analysis_report.html"),
        rreport =os.path.join(config["repo_path"],"{barcode}","analysis_report.html"),
        pcns = os.path.join(config["publish_path"],"{barcode}","consensus_sequences","cns.prompt.txt"),
        rcns = os.path.join(config["repo_path"],"{barcode}","consensus_sequences","cns.prompt.txt"),
        pfkrona = os.path.join(config["publish_path"],"{barcode}","krona.filtered.html"),
        pukrona = os.path.join(config["publish_path"],"{barcode}","krona.unfiltered.html"),
        rfkrona = os.path.join(config["repo_path"],"{barcode}","krona.filtered.html"),
        rukrona = os.path.join(config["repo_path"],"{barcode}","krona.unfiltered.html")
    run:
        shell(
            """
            cp {input.report} {output.preport}
            cp {input.report} {output.rreport}
            cp {input.fkrona} {output.pfkrona}
            cp {input.fkrona} {output.rfkrona}
            cp {input.ukrona} {output.pukrona}
            cp {input.ukrona} {output.rukrona}
            cp {input.cns} {output.pcns}
            cp {input.cns} {output.rcns}
            """)
        try:

            shell("""
            cp {params.cns_path}/*.fasta {params.rcns_path}
            cp {params.cns_path}/*.fasta {params.pcns_path}
            """)
        except:
            print(f"No consensus files to copy for barcode {params.barcode}")


