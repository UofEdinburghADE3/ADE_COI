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
                    if taxid != "0":
                        fw.write(f"{pcent},{subreads},{reads},{rank},{taxid},{name.lstrip()}\n")
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

def get_seq_id_to_tax_map(map_file):
    record_id_to_tax_id = {}
    with open(map_file,"r") as f:
        for l in f:
            l = l.rstrip()
            record,taxid = l.split()
            record_id_to_tax_id[record] = taxid
    return record_id_to_tax_id

def extract_taxid_ref(taxids,map_file,background_file,output_path):
    records = {}

    record_id_to_tax_id = get_seq_id_to_tax_map(map_file)

    for record in SeqIO.parse(background_file, "fasta"):
        try:
            record_taxid = record_id_to_tax_id[record.id]
            if record_taxid in taxids:
                records[record_taxid]= record
        except:
            pass
    valid_taxa = []
    with open(os.path.join(output_path,f"not_animalia.found.txt"),"w") as fw:
        fw.write("taxid\n")
        for record_taxid in taxids:
            if record_taxid not in records:
                fw.write(f"{record_taxid}\n")
            else:
                valid_taxa.append(record_taxid)

    print(f"Number of records found for taxid {record_taxid}: {len(records)}")
    for taxid in records:
        with open(os.path.join(output_path,f"{taxid}.fasta"),"w") as fw:
            record = records[taxid]
            fw.write(f">{record.id}\n{record.seq}\n")
    
    return valid_taxa









