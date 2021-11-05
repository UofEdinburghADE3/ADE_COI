"""
Report structure:

config["report_title"]
config["barcode"]
config["total_reads]
config["total_classified"]
config["total_prop_classified"]

data_for_report["total_classified_svg"]

config["min_length"]
config["max_length"]
config["filtered_reads"]
config["filtered_classified"]
config["filtered_prop_classified"]

Description

Read count: total number of reads for barcode
Reads classified: 
Proportion of reads classified:

Embed
Figure 1: Read length distribution svg

Figure 2: krona plot link

Filtering

Min length (config)
Max length (config)

config["filtered_reads"]
config["filtered_classified"]
config["filtered_prop_classified"]

Read count post filtering: 
Reads classified:
Proportion of reads classified: 

Embed
Figure 3: Read length distribution svg

Figure 4: krona plot link


Result of classification

Minimum number of reads for consensus sequence:
Total number of taxa:

Table 1:
Taxon ID | taxon name | number of reads
xw
"""
import json
import csv
import collections
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date

from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import os
from Bio import SeqIO

from adecoi.utils.log_colours import green,cyan

def count_reads(fastq_file):
    record_dict = SeqIO.index(fastq_file, "fastq")
    return len(record_dict)

def get_read_counts(all_reads, classified_reads, filtered_reads, classified_filtered,data_for_report):
    a_count = count_reads(all_reads)
    ac_count = count_reads(classified_reads)
    f_count = count_reads(filtered_reads)
    fc_count = count_reads(classified_filtered)

    data_for_report["total_reads"] = a_count
    data_for_report["total_classified"] = ac_count
    data_for_report["total_prop_classified"]  = round(a_count/ac_count, 2)
    data_for_report["filtered_reads"] = f_count
    data_for_report["filtered_classified"] = fc_count
    data_for_report["filtered_prop_classified"]  = round(f_count/fc_count, 2)

def load_svgfile(file,key, data_for_report):
    svg = ""
    with open(file,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            svg+=f"{l}\n"
    data_for_report[key] = svg

def load_histograms_svg(hist1, hist2, data_for_report):
    
    load_svgfile(hist1, "histogram1", data_for_report)
    load_svgfile(hist2, "histogram2", data_for_report)

def data_for_table(input_taxa, data_for_report):
    #pcent_reads,sub_reads,reads,rank,taxid,taxon
    data_for_report["taxa_table"] = []
    with open(input_taxa, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report["taxa_table"].append(row)

def make_report(report_to_generate,config,data_for_report,barcode):
    #need to call this multiple times if there are multiple reports wanted
    
    data = json.dumps(data_for_report) 
    print(data)
    template_dir = os.path.abspath(os.path.dirname(config["report_template"]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config["report_template"], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(), 
                    barcode = barcode,
                    data_for_report = data,
                    config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print(green("Generating: ") + f"{report_to_generate}")
        fw.write(buf.getvalue())