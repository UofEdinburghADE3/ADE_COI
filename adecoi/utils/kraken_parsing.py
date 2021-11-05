from Bio import SeqIO
import collections
import os

def parse_report_for_taxa(report,min_reads,out_taxa):
    taxa_dict = {}
    with open(out_taxa,"w") as fw:
        fw.write("pcent_reads,sub_reads,reads,rank,taxid,taxon\n")
        
        with open(report,"r") as f:
            for l in f:
                l = l.rstrip("\n")
                pcent,subreads,reads,rank,taxid,name = l.split("\t")
                if int(reads) > min_reads:
                    fw.write(f"{pcent},{subreads},{reads},{rank},{taxid},{name}\n")
                    taxa_dict[taxid] = int(reads)
    return taxa_dict

def get_reads(taxa_dict,read_file,output_path):
    taxa_read_dict = collections.defaultdict(list)
    for record in SeqIO.parse(read_file,"fastq"):
        try:
            id = record.description.split("taxid|")[1].split()[0]
            if id in taxa_dict:
                taxa_read_dict[id].append(record)
        except:
            pass
    taxa = []
    for taxon in taxa_read_dict:
        reads_found = len(taxa_read_dict[taxon])
        print(f"Found {reads_found} reads in classified file, from reported {taxa_dict[taxon]} reads.")
        if reads_found:
            taxa.append(taxon)
            with open(os.path.join(output_path,f"{taxon}.fastq"),"w") as fw:
                SeqIO.write(taxa_read_dict[taxon],fw,"fastq")
    return taxa

def extract_taxid_ref(taxids,background_file,output_path):
    records = {}
    for record in SeqIO.parse(background_file, "fasta"):
        #5970270|Insecta-Diptera-Chironomidae-||KT605684|kraken:taxid|1721978|gi|930166331
        record_taxid = record.id.split("|")[5]
        if record_taxid in taxids:
            records[record_taxid]= record
    print(f"Number of records found for taxid {record_taxid}: {len(records)}")
    for taxid in records:
        with open(os.path.join(output_path,f"{taxid}.fasta"),"w") as fw:
            record = records[taxid]
            fw.write(f">{record.id}\n{record.seq}\n")








