import os
from Bio import SeqIO
import gzip
import collections
import yaml
import csv
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date

from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO

config = {}
config["index_template"] = "/Users/s1680070/repositories/ADE_COI/adecoi/data/index_template.mako"
config["year"] = "2024"
config["barcodes_map"] = "/Users/s1680070/repositories/ADE_COI/ADE3_2024/barcodes.csv"
config["output_path"] = "/Users/s1680070/repositories/ADE_COI/ADE3_2024"
config["controls"] = []

def make_report(report_to_generate,config,data_for_report,year):
    #need to call this multiple times if there are multiple reports wanted
    
    template_dir = os.path.abspath(os.path.dirname(config["index_template"]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works
    mytemplate = Template(filename=config["index_template"], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(), 
                    year = year,
                    data_for_report = data_for_report,
                    config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, xline %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print("Generating: " + f"{report_to_generate}")
        fw.write(buf.getvalue())


output = os.path.join(config["output_path"],"index.html")

barcode_map = {}
with open(config["barcodes_map"],"r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if  row["barcode"]:
            barcode_map[row["barcode"]]=row["name"].replace("\xa0","")
            
            

data_for_report = []
barcodes_mapped = []
for i in sorted(barcode_map,key= lambda x : int(x)):
    
    if len(i) == 1:
        new_i = f"0{i}"
    else:
        new_i = i
    # print(new_i, barcode_map[i])
    if not os.path.exists(os.path.join(config["output_path"],f"barcode{new_i}")):
        print("This barcode has no report",new_i)
    barcodes_mapped.append(new_i)
    data_for_report.append([barcode_map[i],new_i])

unknown_count = 0
for r,d,f in os.walk(config["output_path"]):
    for dn in sorted(d):
        if "barcode" in dn:
            bc = dn.lstrip("barcode")
            if bc not in barcodes_mapped:
                
                if bc in config["controls"]:
                    data_for_report.append([f"control",bc])
                else:
                    unknown_count+=1
                    data_for_report.append([f"unknown_{unknown_count}",bc])

make_report(output,config,data_for_report,config["year"])




