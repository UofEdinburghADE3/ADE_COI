#!/usr/bin/env python3

from Bio import Entrez

Entrez.email = "aine.otoole@ed.ac.uk"  
handle = Entrez.esearch(db='nucleotide', retmax="40000000", term = ["BOLD[All Fields] AND COI[All Fields]"],idtype="acc")
record = Entrez.read(handle)
count = int(record["Count"])
handle.close()
id_list = record["IdList"]
print(len(id_list), count)
search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))

webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]


batch_size = 100000
out_handle = open("BOLDcoi.fasta", "w")
for start in range(0, count, batch_size):
    end = min(count, start + batch_size)
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
    )
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()