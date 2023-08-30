# YASARA PLUGIN
# TOPIC:       Protein Alignment
# TITLE:       FoldSeek
# AUTHOR:      M.J.L.J. Fürst
# LICENSE:     GPL (www.gnu.org)
# DESCRIPTION: Find Proteins from online databases with foldseek and align
#              Requires foldseek to be installed.

"""
MainMenu: Analyze
  PullDownMenu: Align
    Submenu after Molecules with FatCat: FoldSeek
      CustomWindow: Foldseek parameters
        Width: 600
        Height: 380
        Text:        X= 20,Y=48,Text="Choose database:"
        RadioButtons:Options=3,Default=1
                     X=20,Y=68,Text="PDB"
                     X=20,Y=103,Text="AlphaFold Swissprot"
                     X=20,Y=138,Text="AlphaFold UniProt50 (slow)" 
        CheckBox:    X= 265,Y=48,Text="Delete non-homologous chains",Default=Yes
        NumberInput: X= 265,Y= 90,Text="Number of structures to retrieve",Default=20,Min=1,Max=1000
        CheckBox:    X= 390,Y=107,Text="All",Default=No
        TextInput:   X=20, Y= 190, Text='Output folder (current folder if empty)', 550, 100
        TextInput:   X=20, Y= 260, Text='Custom option flags (check manual)', 550, 100
        Button:      X=281,Y=330,Text="_O_ K"
      ObjectSelectionWindow: Select objects to find tunnels for
      Request: foldseek
"""

### custom functions
def check_and_install_module(module_name):
    try:
        importlib.import_module(module_name)
        return True
    except ImportError:
        return False
    

def stop_plg(message, fs=True):
    PrintCon()
    if fs:
        Print(f'Ran FoldSeek plugin on object {selection[1].objects} searching the {["PDB", "AF SwissProt", "AF Uniprot"][selection[0].radiobutton -1]} database for the {get} closest structural homologs.\nTime for execution of script: {str(timedelta(seconds=float("{:.1f}".format(time.perf_counter()  - start_time))))}')
    ShowMessage(message)
    Wait('Continuebutton')
    HideMessage()
    plugin.end()

def find_subfolder(target_folder_name, parent_folder_name, start_directory):
    for dirpath, dirnames, filenames in os.walk(start_directory):
        if parent_folder_name in dirnames:
            parent_folder_path = os.path.join(dirpath, parent_folder_name)
            for sub_dirpath, sub_dirnames, sub_filenames in os.walk(parent_folder_path):
                if target_folder_name in sub_dirnames:
                    return os.path.join(sub_dirpath, target_folder_name)
    return None


def convert_cif_to_pdb(cif_file_path, pdb_file_path):
    # Create a CIF parser and parse the CIF file
    cif_parser = MMCIFParser()
    structure = cif_parser.get_structure("structure", cif_file_path)

    # Create a PDBIO object to write the PDB file
    pdb_io = PDBIO()

    # Save the structure as a PDB file
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_file_path)


def download_pdb(server, id, outfile):
    url = f'{server}/{id}.pdb'
    response = requests.get(url)
    if response.status_code == 200:
        with open(outfile, 'wb') as file:
            file.write(response.content)
    elif response.status_code == 404:
        url = f'{server}/{id}.cif'
        response = requests.get(url)
        if response.status_code == 200:
            with open(outfile.replace('pdb', 'cif'), 'wb') as file:
                file.write(response.content)
            convert_cif_to_pdb(outfile.replace('pdb', 'cif'), outfile)
        else:
            print(f"Status code {response.status_code} caused failure to download file {url} ")
    else:
        print(f"Status code {response.status_code} caused failure to download file {url} ")

from yasara import *

if request == 'foldseek':
    ### dependency checks
    import importlib
    import subprocess
    import sys

    print(f'python interpreter is: {sys.executable}')
    required_modules = {'pandas': 'pandas', # Import name before colon, pip installation name after colon
                        'requests': 'requests',
                        'gzip': 'gzip',
                        'Bio': 'Bio'} 
    missing_modules = [module for module in required_modules if not check_and_install_module(module)]
    if missing_modules:
        install_choice =\
            ShowWin("Custom","Missing software",600,180,
                    "Text",         20, 48,"The following python modules weren't found on your system:",
                    "Text",         20, 73, ",".join(missing_modules),
                    "Text",         20, 98,"Do you want to install these automatically?",
                    "Button",      250,130,"No",
                    "Button",      350,130,"Yes")[0]
        
        if install_choice == 'Yes':
            python_executable = sys.executable
            for module in missing_modules:
                install_name = required_modules[module]
                subprocess.call([python_executable, '-m', 'pip', 'install', install_name])
        else:
            ShowMessage("Some python modules are missing, try installing manually. Exiting.")
            plugin.end()

    import os
    import pandas as pd
    import requests
    import gzip
    from collections import OrderedDict
    from Bio.PDB import MMCIFParser, PDBIO
    import shutil
    import re
    from datetime import timedelta

    ### process user parameters
    # which database of foldseek to use and where is it? where can files be downloaded?
    if selection[0].radiobutton == 1:
        fsdb = '/usr/local/bin/foldseek/pd'
        dlserver = 'https://files.rcsb.org/view'
    elif selection[0].radiobutton == 2:
        fsdb = '/usr/local/bin/foldseek/sp'
        dlserver = 'https://alphafold.ebi.ac.uk/files'
    elif selection[0].radiobutton == 3:
        fsdb = '/usr/local/bin/foldseek/up50'
        dlserver = 'https://alphafold.ebi.ac.uk/files'

    # where to save the output
    if os.path.exists(selection[0].text[0]):
        outputdir = selection[0].text[0]
    else:
        if selection[0].text[0] == '':
            outputdir = PWD()
        else:
            stop_plg('Error: output directory doesn\'t exist.')

    # how many structures to get
    if selection[0].checkbox[1]:
        confirm = ShowWin('Custom', 'Warning', 400, 170, 
                        'Text', 20, 48, "Downloading all homologs can take a long",
                        'Text', 20, 78, "time and cause memory overload. Continue?",
                        'Button', 130, 120, "Yes",
                        'Button', 250, 120, "No")
        if confirm[0] == 'Yes':
            get = 'all'
        else:
            stop_plg('Aborted.')
    else:
        get = selection[0].number[0]


    ### get foldseek executable by trying `which`, then chcecking conda bin dir, then searching everywhere
    fs = shutil.which('foldseek')
    if not fs:
        if os.path.exists(os.path.dirname(os.environ.get('CONDA_EXE', None)) + os.path.sep + 'foldseek'):
            fs = os.path.dirname(os.environ.get('CONDA_EXE', None)) + os.path.sep + 'foldseek'
        else:
            fs = find_subfolder('foldseek', 'bin', '/')
            if fs:
                fs = fs + '/bin/foldseek'

    if fs:
        version = re.search('(?<=Version: ).*', subprocess.run(f'{fs} -h', shell=True, capture_output=True, text=True).stdout).group()
        print(f'Using foldseek version {version} at {fs}')
    else:
        stop_plg('Error: Foldseek not found. Aborting.')

    # if installed version is >=6, a more memory effcient option is available
    # if int(version[0]) >= 6:
    if version[0] == 'b' or int(version[0]) >= 6:
        options = '--prefilter-mode 1 ' + selection[0].text[1]
    else:
        options = selection[0].text[1]

    ### Begin of actual script
    start_time = time.perf_counter()
    Console('OFF')

    # make output dirs if necessary
    dirs = [outputdir + os.path.sep + 'tmp', outputdir + os.path.sep + 'q', outputdir + os.path.sep + 'alns', outputdir + os.path.sep + 'hits']
    for dir in dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

    target = [selection[1].object[j].number.inyas for j in range(selection[1].objects)][0]
    target_name = NameObj(target)[0]

    # create the query pdb and run foldseek
    SavePDB(target, outputdir + os.path.sep + 'q' + os.path.sep + target_name + '_fsquery.pdb')
    ShowMessage(f'Running FoldSeek on object {target}, please wait.')
    fs_command = f'{fs} easy-search {outputdir}/q/{target_name}_fsquery.pdb {fsdb} {outputdir}/alns/{target_name}_aln {outputdir}/tmp --remove-tmp-files true {options}'
    Print('Command:\n' + fs_command)
    Wait(1)
    noerr = subprocess.run(fs_command, 
                shell=True, capture_output=True, text=True)
    for line in noerr.stdout.split('\n'):
        print(line)

    # ensure it has worked
    if noerr.returncode > 0 or not os.path.exists(f'{outputdir}/alns/{target_name}_aln'):
        stop_plg('Error: Something went wrong, check the console.')

    # load hits from file, get unique pdbs and number of requested results
    try:
        hits = pd.read_table(f'{outputdir}/alns/{target_name}_aln', header=None)[1]
    except pd.errors.EmptyDataError:
        stop_plg('No homologs were found!')

    # check if a header was present
    if hits.iloc[0] == 'target':
        hits = hits.iloc[1:].tolist()  # Discard the first row
    else:
        hits = hits.tolist()
    
    hits = list(OrderedDict.fromkeys(hits))
    if get != 'all':
        hits = hits[:get]

    ShowMessage(f'Getting {len(hits)} structures, please wait.')
    Wait(1)

    ### download pdbs
    # write hit file
    if dlserver == 'https://files.rcsb.org/view':
        hit_pdbs = [x.split('_')[0] for x in hits]
        hit_mols = [x.split('_')[1] for x in hits]
    elif dlserver == 'https://alphafold.ebi.ac.uk/files':
        hit_pdbs = hits
        hit_mols = ['A' for x in range(len(hits))]
        
    with open(f'{outputdir}/hits/{target_name}_hits', 'w') as file:
        file.write('\n')
        for i in range(len(hit_pdbs)):
            file.write(f" {hit_pdbs[i]}   {hit_mols[i]}\n")

    # download with rsync if more than 20, otherwise and non-rsynced manually
    if get == 'all' or get > 20:
        with open(f'{outputdir}/hits/{target_name}_hits_rslist', 'w') as file:
            for i in range(len(hit_pdbs)):
                file.write(f"pdb{hit_pdbs[i]}.ent.gz\n")   

        noerr = subprocess.run(f'rsync -rlptL -v -z --delete --relative --files-from={outputdir}/hits/{target_name}_hits_rslist --port=33444 rsync.wwpdb.org::ftp/data/structures/all/pdb/ {outputdir}/hits/', 
                shell=True, capture_output=True, text=True)
        for line in noerr.stdout.split('\n'):
            print(line)
        
        for i in range(len(hit_pdbs)):
            if os.path.exists(f'{outputdir}/hits/pdb{hit_pdbs[i]}.ent.gz'):
                print(f'file {hit_pdbs[i]} was rsynced, extracting.')
                with gzip.open(f'{outputdir}/hits/pdb{hit_pdbs[i]}.ent.gz', 'rb') as input_file:
                    with open(f'{outputdir}/hits/{i+1}_{hit_pdbs[i]}.pdb', 'wb') as output_file:
                        output_file.write(input_file.read())
                os.remove(f'{outputdir}/hits/pdb{hit_pdbs[i]}.ent.gz')
            else:
                print(f'file {hit_pdbs[i]} wasn\'t rsynced, getting it manually.')
                download_pdb(dlserver, hit_pdbs[i], f'{outputdir}/hits/{i+1}_{hit_pdbs[i]}.pdb')
    else:
        for i in range(get):
            download_pdb(dlserver, hit_pdbs[i], f'{outputdir}/hits/{i+1}_{hit_pdbs[i]}.pdb')
        

    for i in range(len(hit_pdbs)):
        ShowMessage(f'Loading homolog {i + 1} / {len(hit_pdbs)}')
        Wait(1)
        new = LoadPDB(f'{outputdir}/hits/{i+1}_{hit_pdbs[i]}.pdb')
        # could do for all structures in nmr bundles
        # [HideObj(i) for i in new]
        # [DelMol(f'obj {i} mol !{hit_mols[i]}') for i in new]
        # [AlignObj(i, target, 'sheba') for i in new]
        [DelObj(i) for i in new if i > min(new)]
        HideObj(min(new))
        if selection[0].checkbox[0]:
            DelMol(f'obj {min(new)} mol !{hit_mols[i]}')
        AlignObj(min(new), target, 'sheba')

    if dlserver == 'https://files.rcsb.org/view':
        img = MakeImage("Buttons",topcol="None",bottomcol="None")
        ShowImage(img,alpha=66,priority=0)
        PrintImage(img)

        def ShowButtons():
            Font("Arial",height=14,color="black")

            ShowButton("Show info",x='12%', y='65%',color="White", height=40)
            ShowButton("Open PDB site",x='12%', y='70%',color="White", height=40)
            ShowButton("Open primary article",x='12%', y='75%',color="White", height=40)
            ShowButton("Exit",x='12%', y='80%',color="Red", height=40)


        ShowButtons()
    stop_plg('Finished.')

elif request == 'Exit':
    SaveSce('ExitFoldSeek.sce')
    Clear()
    LoadSce('ExitFoldSeek.sce')

elif request == 'Showinfo':
    import os
    import re
    Console('OFF')

    def find_matching_files(folder_path, input_strings):
        matching_files = []
        for file_name in os.listdir(folder_path):
            for input_str in input_strings:
                if file_name.startswith(input_str) and file_name.endswith('_aln'):
                    matching_files.append(os.path.join(folder_path, file_name))
        return matching_files
    
    aln_file = find_matching_files(f'{PWD()}/alns/', NameObj('ALL'))[0]
    target = ListObj(ShowWin('ObjectSelection', 'choose PDB')[0], format='OBJNAME')[0]
    target = re.sub('^[0-9]+_', '', target)
    vals = [line for line in open(f'{PWD()}/alns/5m10_aln') if target in line][0].split('\t')
    cols = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'.split(',')
    
    print(f"Info from {aln_file}:")
    for token, value in zip(cols, vals):
        print(f"{token.ljust(max(len(token) for token in cols))}: {value}")

elif request == 'OpenPDBsite':
    import re
    Console('off')

    target = ListObj(ShowWin('ObjectSelection', 'choose PDB')[0], format='OBJNAME')[0]
    target = re.sub('^[0-9]+_', '', target)
    ShowURL(f'https://www.rcsb.org/structure/{target}')

elif request == 'Openprimaryarticle':
    import re
    import requests
    Console('off')

    target = ListObj(ShowWin('ObjectSelection', 'choose PDB')[0], format='OBJNAME')[0]
    target = re.sub('^[0-9]+_', '', target)
    url = f'https://www.rcsb.org/structure/{target}'

    response = requests.get(url)

    if response.status_code == 200:
        doi = re.findall('https:\/\/doi\.org\/[\w\/.-]+', response.text)
        doi = [x for x in list(set(doi)) if not x.endswith('pdb')]
        print(f"found {len(doi)} article{'s' * (len(doi) != 1)} for this structures:")
        for d in doi:
            Print(d)
        if len(doi) == 0:
            stop_plg(f"Could not identify doi from this pdb site: {url}", fs=False)
        else:
            ShowURL(doi[0])
    else:
        stop_plg(f"Request failed with status code: {response.status_code}", fs=False)
        plugin.end()

plugin.end()
