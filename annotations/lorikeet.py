import sys
import os
from distutils.dir_util import copy_tree
import shutil
import numpy as np
import subprocess, linecache

htmlpage = '''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
 <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    
    <title>Lorikeet Spectrum Viewer</title>
    
    <!--[if IE]><script language="javascript" type="text/javascript" src="js/excanvas.min.js"></script><![endif]-->
    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script>
    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.4/jquery-ui.min.js"></script>
    <script type="text/javascript" src="js/jquery.flot.js"></script>
    <script type="text/javascript" src="js/jquery.flot.selection.js"></script>
    
    <script type="text/javascript" src="js/specview.js"></script>
    <script type="text/javascript" src="js/peptide.js"></script>
    <script type="text/javascript" src="js/aminoacid.js"></script>
    <script type="text/javascript" src="js/ion.js"></script>
    
    <link REL="stylesheet" TYPE="text/css" HREF="css/lorikeet.css">
    
</head>

<body>

<!-- PLACE HOLDER DIV FOR THE SPECTRUM -->
<div id="lorikeet1"></div>

<script type="text/javascript">

$(document).ready(function () {

	/* render the spectrum with the given options */
	$("#lorikeet1").specview({sequence: sequence, 
								charge: charge,
								massError: 0.02,
								precursorMz: precursorMz,
								variableMods: varMods, 
								//ctermMod: ctermMod,
								peaks: peaks
								});	
});

'''

# Parsing an MS2 spectrum (title) from an MGF file (mgf)
# As spectrum files can be quite large, I use findstr to read the lines first
# This only really speeds up things if the next query is on the same MGF file
# TODO: We should think about how to optimize this for ionbot.cloud
def get_spectrum(mgf, scan):
    line = subprocess.check_output(['findstr', '/N', "SCANS=%i"%scan, mgf])

    line = int(line.decode("utf-8").split(":")[0])
    spectrum = "["
    while True:
        c = linecache.getline(mgf, line)
        c = c.rstrip()
        if c == "": 
            line+=1
            continue
        if "END IONS" in c: break
        if "PEPMASS=" in c:
            parent_mz = c[8:]
        if "CHARGE" in c:
            charge = c[7:9].replace("+","")
        if not "=" in c:
            tmp = c.split(" ")
            spectrum += "[%s,%s],"%(tmp[0],tmp[1])
        line+=1
    spectrum = spectrum[:-1]
    spectrum += "]"
    return spectrum, charge, parent_mz

# Here the "matched_peptide" and "modifications" columns
# in the ionbot result file are passed to create the data
# for the varMods javascript variable
def get_varmods(peptide, modifications):
    mods = []
    tmp = modifications.split("|")
    for i in range(0,len(tmp),2):
        mod_pos = int(tmp[i])
        mods.append("{index: %i, modMass: %s, aminoAcid: '%s'}"%(mod_pos,tmp[i+1],peptide[mod_pos-1]))
    return mods

def generate_html(annotations_folder,mgf_folder,mgf_file,scan,df):
    tmp = df[(df["spectrum_file"]==mgf_file)&(df["scan"]==scan)]
    sequence = tmp["matched_peptide"].values[0]
    modifications = tmp["modifications_delta"].values[0]
    spectrum, charge, parent_mz = get_spectrum(mgf_folder+mgf_file, scan)
    if spectrum == "]":
        print("spectrum not found")

    varmods_list = []
    if modifications != "0|":
        varmods_list = get_varmods(sequence, modifications)

    if os.path.exists(annotations_folder) == False:
        os.mkdir(annotations_folder)
        shutil.copytree("css",os.path.join(annotations_folder,"css"))
        shutil.copytree("js",os.path.join(annotations_folder,"js"))
    
        
    html_filename = annotations_folder + sequence+"_"+str(scan)+".html"
    with open(html_filename,'w') as f:
        f.write(htmlpage+'\n')
        f.write('var sequence = "%s";\n'%sequence)
        f.write('var peaks = %s;\n'%spectrum)
        f.write('var charge = %s;\n'%charge)
        #f.write('var precursorMz = %s;\n'%parent_mz)
        f.write('var precursorMz = 0;\n')
        f.write('var varMods = [];\n')
        for i,mod in enumerate(varmods_list):
            f.write("varMods[%i] = %s\n"%(i,mod))
        f.write('</script></body></html>\n')
    return html_filename