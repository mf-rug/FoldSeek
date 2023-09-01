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
    

def stop_plg(message, fs=True, start_time=None):
    if start_time == None:
        start_time = time.perf_counter()
    PrintCon()
    if fs:
        Print(f'Ran FoldSeek plugin on object {selection[0].objects} searching the {databases[int(db) -1]} database for the {n_get} closest structural homologs.\nTime for execution of script: {timedelta(seconds=float("{:.1f}".format(time.perf_counter() - start_time)))}')
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


def get_pdbs(pdb_ids, dlserver, dir, target):
    # for pdb: download with rsync if more than 20, otherwise and non-rsynced manually
    if dlserver == 'https://files.rcsb.org/view' and len(pdb_ids) > 20:
        with open(os.path.join(dir, f'{target}_hits_rslist'), 'w') as file:
            for i in range(len(pdb_ids)):
                file.write(f"pdb{pdb_ids[i]}.ent.gz\n")   

        noerr = subprocess.run(f'rsync -rlptL -v -z --delete --relative --files-from={os.path.join(dir, f"{target}_hits_rslist")} --port=33444 rsync.wwpdb.org::ftp/data/structures/all/pdb/ {outputdir}/hits/', 
                shell=True, capture_output=True, text=True)
        for line in noerr.stdout.split('\n'):
            print(line)
        
        for i in range(len(pdb_ids)):
            if os.path.exists(f'{dir}/pdb{pdb_ids[i]}.ent.gz'):
                print(f'file {pdb_ids[i]} was rsynced, extracting.')
                with gzip.open(f'{dir}/pdb{pdb_ids[i]}.ent.gz', 'rb') as input_file:
                    with open(f'{dir}/{i+1}_{pdb_ids[i]}.pdb', 'wb') as output_file:
                        output_file.write(input_file.read())
                os.remove(f'{dir}/pdb{pdb_ids[i]}.ent.gz')
            else:
                print(f'file {pdb_ids[i]} wasn\'t rsynced, getting it manually.')
                download_pdb(dlserver, pdb_ids[i], f'{dir}/{i+1}_{pdb_ids[i]}.pdb')
    else:
        for i in range(n_get):
            download_pdb(dlserver, pdb_ids[i], f'{dir}/{i+1}_{pdb_ids[i]}.pdb')


def download_pdb(server, id, outfile):
    url = f'{server}/{id}.pdb'
    response = get(url)
    if response.status_code == 200:
        with open(outfile, 'wb') as file:
            file.write(response.content)
    elif response.status_code == 404:
        url = f'{server}/{id}.cif'
        response = get(url)
        if response.status_code == 200:
            with open(outfile.replace('pdb', 'cif'), 'wb') as file:
                file.write(response.content)
            convert_cif_to_pdb(outfile.replace('pdb', 'cif'), outfile)
        else:
            print(f"Status code {response.status_code} caused failure to download file {url} ")
    else:
        print(f"Status code {response.status_code} caused failure to download file {url} ")


def load_pdbs(pdbs, mols, outdir, del_homologs):
    for i in range(len(pdbs)):
        ShowMessage(f'Loading homolog {i + 1} / {len(pdbs)}')
        Wait(1)
        new = LoadPDB(f'{outdir}{i+1}_{pdbs[i]}.pdb')
        # could do for all structures in nmr bundles
        # [HideObj(i) for i in new]
        # [DelMol(f'obj {i} mol !{mols[i]}') for i in new]
        # [AlignObj(i, target, 'sheba') for i in new]
        [DelObj(i) for i in new if i > min(new)]
        HideObj(min(new))
        if del_homologs:
            DelMol(f'obj {min(new)} mol !{mols[i]}')
        AlignObj(min(new), target, 'sheba')


def ShowButtons(dlserver):
    img = MakeImage("Buttons",topcol="None",bottomcol="None")
    ShowImage(img,alpha=66,priority=0)
    PrintImage(img)
    Font("Arial",height=14,color="black")

    ShowButton("Show info",x='12%', y='65%',color="White", height=40)
    ShowButton("Open structure entry",x='12%', y='70%',color="White", height=40)
    if dlserver == 'https://files.rcsb.org/view':
        ShowButton("Open primary article",x='12%', y='75%',color="White", height=40)
        ShowButton("Exit",x='12%', y='80%',color="Red", height=40)
    else
        ShowButton("Exit",x='12%', y='75%',color="Red", height=40)


from yasara import *
import os
import re
from requests import get, post
from tempfile import gettempdir

Console('OFF')

if request == 'foldseek':
    ### dependency checks
    import importlib
    import subprocess
    import sys
    import shutil
    from datetime import timedelta

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
        online_only = False


    ### ONLINE FOLDSEEEK
    if online_only or os.path.exists(os.path.join(gettempdir(), 'fs_plg_runonline')):
        # this is a bit of an ugly solution to run foldseek online if the user requests it.
        if os.path.exists(os.path.join(gettempdir(), 'fs_plg_runonline')):
            Print('User specified online version this time!')
            os.remove(os.path.join(gettempdir(), 'fs_plg_runonline'))

        start_time = time.perf_counter()
        db, del_homologs, n_get, get_all, out_dir =\
            ShowWin("Custom", "Foldseek Webserver Parameters", 600, 355,
                "Text", 20, 48, "Choose database:",
                "RadioButtons", 5, 1,
                20,68,"PDB",
                20,103,"AlphaFold Swissprot",
                20,138,'AlphaFold UniProt50',
                20,173,'AlphaFold Proteome',
                20,208,'MGnify/ESM30',
                "CheckBox", 265, 48, "Delete non-homologous chains", True,
                "NumberInput", 265, 90, "Number of structures to retrieve", 20, 1, 1000,
                "CheckBox", 390, 107, "All", False,
                "TextInput", 20, 255, "Output folder (current folder if empty)", 550, 100,
                "Button", 281, 315, "_O_ K")
        
        databases = ['afdb50', 'afdb-swissprot', 'afdb-proteome', 'cath50', 'mgnify_esm30', 'pdb100', 'gmgcl_id']
        database = databases[db -1]
        print(f'Selected database {database}')

        from time import sleep
        import sys
        import tarfile
        import datetime
        import pandas as pd

        # where to save the output
        if os.path.exists(out_dir):
            outputdir = out_dir
        else:
            if out_dir == '':
                outputdir = PWD()
            else:
                stop_plg('Error: output directory doesn\'t exist.', start_time=start_time)

        # how many structures to get
        if get_all:
            confirm = ShowWin('Custom', 'Warning', 400, 170, 
                            'Text', 20, 48, "Downloading all homologs can take a long",
                            'Text', 20, 78, "time and cause memory overload. Continue?",
                            'Button', 130, 120, "Yes",
                            'Button', 250, 120, "No")
            if confirm[0] == 'Yes':
                n_get = 'all'
            else:
                stop_plg('Aborted.', start_time=start_time)

        target = [selection[0].object[j].number.inyas for j in range(selection[0].objects)][0]
        target_name = NameObj(target)[0]

        # make output dirs if necessary
        target_out_dir = os.path.join(outputdir, 'fs_web_hits', target_name)
        dirs = [os.path.join(outputdir, 'fs_web_hits'), target_out_dir]
        for dir in dirs:
            if not os.path.exists(dir):
                os.mkdir(dir)

        # Get the current date and time
        datetime_str = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        fquery = os.path.join(target_out_dir, target_name + '_' + datetime_str + '_fsquery.pdb')
        SavePDB(target, fquery)

        # submit a new job
        ShowMessage('Submitting job to FoldSeek webserver.')
        Wait(1)
        with open(fquery, 'r') as pdb_file:
            pdb_content = pdb_file.read()
            ticket = post('https://search.foldseek.com/api/ticket', 
                        files={'q': (pdb_content, pdb_content, 'application/octet-stream')},
                        data={
                            'mode' : '3diaa',
                            'database[]' : databases,
                        }).json()

        # poll until the job was successful or failed
        repeat = True
        wait_time = time.perf_counter()
        while repeat:
            ShowMessage(f'Waiting for result from FoldSeek server ({round(time.perf_counter()  - wait_time,0)}s)')
            Wait(1)
            status = get('https://search.foldseek.com/api/ticket/' + ticket['id']).json()
            if status['status'] == "ERROR":
                # handle error
                sys.exit(0)

            # wait a short time between poll requests
            sleep(10)
            repeat = status['status'] != "COMPLETE"
        
        ShowMessage(f'Receiving results ({round(time.perf_counter()  - wait_time,0)}s)')
        Wait(1)
        # get results in JSON format
        result = get('https://search.foldseek.com/api/result/' + ticket['id'] + '/0').json()

        # download result to file
        download = get('https://search.foldseek.com/api/result/download/' + ticket['id'], stream=True)
        web_out_file = os.path.join(target_out_dir, 'result.tar.gz')
        
        with open(web_out_file, 'wb') as fd:
            for chunk in download.iter_content(chunk_size=128):
                fd.write(chunk)

        if os.path.exists(web_out_file):
            file = tarfile.open(web_out_file)
            file.extractall(target_out_dir)            
            file.close()
            out_files = [x for x in os.listdir(target_out_dir) if x.startswith('ali') and x.endswith('m8')]
            for file in out_files:
                os.rename(os.path.join(target_out_dir, file), os.path.join(target_out_dir, target_name + '_' + datetime_str + '_' + file))
        else:
            stop_plg(f"An error occured while running obj {target} the FoldSeek server.", start_time=start_time)

        alns = os.path.join(target_out_dir, target_name + '_' + datetime_str + '_' + [x for x in out_files if database in x][0])

        server_output_column_names = ["query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","prob","evalue","bits","qlen","tlen","qaln","taln","tca","tseq","taxid","taxname"]
        hits_df = pd.read_table(alns, header=None, names=server_output_column_names)

        hits = hits_df['target'].tolist()
        hits = list(set([x.split(' ')[0] for x in hits]))

        if n_get != 'all':
            hits = hits[:n_get]

        ShowMessage(f'Getting {len(hits)} structures, please wait.')
        Wait(1)

        ### download pdbs
        # write hit file
        if 'afdb' in database:
            dlserver = 'https://alphafold.ebi.ac.uk/files'
            hit_pdbs = hits
            hit_mols = ['A' for x in range(len(hits))]
        elif database == 'pdb100': 
            dlserver = 'https://files.rcsb.org/view'
            hit_pdbs = [x.split('_')[0] for x in hits]
            hit_mols = [x.split('_')[1] for x in hits]
        elif database == 'mgnify_esm30':
            dlserver = 'https://api.esmatlas.com/fetchPredictedStructure'
            hit_pdbs = [x.split('.')[0] for x in hits]
            hit_mols = ['A' for x in range(len(hits))]
        
        if not os.path.exists(f'{target_out_dir}/hits/'):
            os.mkdir(f'{target_out_dir}/hits/')
        
        with open(os.path.join(target_out_dir, 'hits', target_name + '_' + datetime_str + '_' + 'hits'), 'w') as file:
            file.write('\n')
            for i in range(len(hit_pdbs)):
                file.write(f" {hit_pdbs[i]}   {hit_mols[i]}\n")
             
        get_pdbs(hit_pdbs, dlserver=dlserver, dir=f'{target_out_dir}/hits/', target=target_name)
        if dlserver == 'https://api.esmatlas.com/fetchPredictedStructure':
            import fileinput
            pattern = r'(^ATOM +[0-9]+)( *[A-Z0-9]+) (.*)'
            replacement = r'\1 \2\3'
            for i, hit in enumerate(hit_pdbs):
                # Construct the full path to the file
                file_path = os.path.join(os.path.join(target_out_dir, 'hits'), f'{i + 1}_{hit}.pdb')

                # check if pbd file has a formatting issue
                chk_pattern = '^ATOM + [0-9]+ +(?=[^ ] )'
                chk_line = next((re.search(chk_pattern, line).group() for line in open(file_path) if re.search(chk_pattern, line)), None)
                if len(''.join(reversed(chk_line)).lstrip()) - len(chk_line) < 0:
                    print(f'Found formatting issue in file {file_path}; auto-correcting.')
                    # Process the file in-place
                    with fileinput.FileInput(file_path, inplace=True, backup='.bak') as file:
                        for line in file:
                            modified_line = re.sub(pattern, replacement, line)
                            print(modified_line, end='')

        load_pdbs(hit_pdbs, hit_mols, f'{target_out_dir}/hits/', del_homologs= del_homologs)
  
        ShowMessage('Done.')

    ### LOCAL / OFFLINE FOLDSEEEK
    else:
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
                    ShowMessage(f'Installing {install_name}, please wait.')
                    cmd = [python_executable, '-m', 'pip', 'install', install_name]
                    no_err = subprocess.run(subprocess.list2cmdline(cmd),shell=True, capture_output=True, text=True)
                    for line in no_err.stdout.split('\n'):
                        print(line)
                    if no_err.stderr != '':
                        print(no_err.stderr)
                missing_modules = [module for module in required_modules if not check_and_install_module(module)]
                if missing_modules:
                    stop_plg("Some python modules couldn't be installed, see the console. Try installing manually. Exiting.")
                else:
                    ShowMessage('Installation successful.')
            else:
                ShowMessage("Some python modules are missing, try installing manually. Exiting.")
                plugin.end()

        import pandas as pd
        import gzip
        from collections import OrderedDict
        from Bio.PDB import MMCIFParser, PDBIO

        button_click, db, del_homologs, n_get, get_all, out_dir, flags  =\
            ShowWin("Custom", "Local FoldSeek Parameters", 600, 380,
                "Text", 20, 48, "Choose database:",
                "RadioButtons", 3, 1,
                20,68,"PDB",
                20,103,"AlphaFold Swissprot",
                20,138,'AlphaFold UniProt50 (slow)',
                "CheckBox", 265, 48, "Delete non-homologous chains", True,
                "NumberInput", 265, 90, "Number of structures to retrieve", 20, 1, 1000,
                "CheckBox", 390, 107, "All", False,
                "TextInput", 20, 190, "Output folder (current folder if empty)", 550, 100,
                "TextInput", 20, 260, "Custom option flags (check manual)", 550, 100,
                "Button", 350, 330, "Switch to Online FoldSeek",
                "Button", 150, 330, "OK")

        databases = ["PDB", "AF SwissProt", "AF Uniprot"]
        if button_click != 'OK':
            with open(os.path.join(gettempdir(), 'fs_plg_runonline'), 'w') as file:
                file.write('run_online')
            stop_plg('Run the plugin again, it will start in online mode this time.')

        ### process user parameters
        # if installed version is >=6, a more memory effcient option is available
        # if int(version[0]) >= 6:
        if version[0] == 'b' or int(version[0]) >= 6:
            options = '--prefilter-mode 1 ' + flags
        else:
            options = flags

        # which database of foldseek to use and where is it? where can files be downloaded?
        if db == 1:
            fsdb = '/usr/local/bin/foldseek/pd'
            dlserver = 'https://files.rcsb.org/view'
        elif db == 2:
            fsdb = '/usr/local/bin/foldseek/sp'
            dlserver = 'https://alphafold.ebi.ac.uk/files'
        elif db == 3:
            fsdb = '/usr/local/bin/foldseek/up50'
            dlserver = 'https://alphafold.ebi.ac.uk/files'

        # where to save the output
        if os.path.exists(out_dir):
            outputdir = out_dir
        else:
            if out_dir == '':
                outputdir = PWD()
            else:
                stop_plg('Error: output directory doesn\'t exist.', start_time)

        # how many structures to get
        if get_all:
            confirm = ShowWin('Custom', 'Warning', 400, 170, 
                            'Text', 20, 48, "Downloading all homologs can take a long",
                            'Text', 20, 78, "time and cause memory overload. Continue?",
                            'Button', 130, 120, "Yes",
                            'Button', 250, 120, "No")
            if confirm[0] == 'Yes':
                n_get = 'all'
            else:
                stop_plg('Aborted.')

        ### Begin of actual script
        start_time = time.perf_counter()

        # make output dirs if necessary
        dirs = [outputdir + os.path.sep + 'tmp', outputdir + os.path.sep + 'q', outputdir + os.path.sep + 'alns', outputdir + os.path.sep + 'hits']
        for dir in dirs:
            if not os.path.exists(dir):
                os.mkdir(dir)

        target = [selection[0].object[j].number.inyas for j in range(selection[0].objects)][0]
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
            stop_plg('Error: Something went wrong, check the console.', start_time=start_time)

        # load hits from file, get unique pdbs and number of requested results
        try:
            hits = pd.read_table(f'{outputdir}/alns/{target_name}_aln', header=None)[1]
        except pd.errors.EmptyDataError:
            stop_plg('No homologs were found!', start_time=start_time)

        # check if a header was present
        if hits.iloc[0] == 'target':
            hits = hits.iloc[1:].tolist()  # Discard the first row
        else:
            hits = hits.tolist()
        
        hits = list(OrderedDict.fromkeys(hits))
        if n_get != 'all':
            hits = hits[:n_get]

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
             
        get_pdbs(hit_pdbs, dlserver=dlserver, dir=f'{outputdir}/hits/', target=target_name)

        load_pdbs(hit_pdbs, hit_mols, f'{outputdir}/hits/', del_homologs= del_homologs)

    ShowButtons(dlserver)
    stop_plg('Finished.', start_time=start_time)

elif request == 'Exit':
    SaveSce('ExitFoldSeek.sce')
    Clear()
    LoadSce('ExitFoldSeek.sce')

elif request == 'Showinfo':

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

elif request == 'Openstructureentry':
    target = ListObj(ShowWin('ObjectSelection', 'choose PDB')[0], format='OBJNAME')[0]
    target = re.sub('^[0-9]+_', '', target)
    if re.search('^MG', target):
        ShowURL(f'https://esmatlas.com/resources/detail/{target}')
    elif re.search('^AF-', target):
        ShowURL(f'https://www.rcsb.org/structure/{target}')
    else:
        ShowURL(f'https://alphafold.ebi.ac.uk/entry/{target}')

elif request == 'Openprimaryarticle':
    target = ListObj(ShowWin('ObjectSelection', 'choose PDB')[0], format='OBJNAME')[0]
    target = re.sub('^[0-9]+_', '', target)
    url = f'https://www.rcsb.org/structure/{target}'

    response = get(url)

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
