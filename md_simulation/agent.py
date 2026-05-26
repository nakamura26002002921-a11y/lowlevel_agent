import sys,os,requests,glob
import subprocess
import re
import shlex
import shutil
import json
from pathlib import Path
from pymol import cmd
from contextlib import contextmanager
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

BASE_DIR = os.path.dirname(__file__)
SHELL_DIR = str(Path(__file__).parent.parent / "shellcommands")

def initialization(base_path):
    if not base_path:
        raise ValueError("PATH is required")
    playbook_path = os.path.join(BASE_DIR, "create_dirs.yml")
    cmd = ["ansible-playbook", playbook_path, "-i", "localhost,", "-c", "local", "--extra-vars", f"PATH={shlex.quote(str(base_path))}"]
    subprocess.run(cmd, check=True)

def copy_mdp(base_path, reference_mdp_path):
    if not base_path:
        raise ValueError("PATH is required.")
    if not reference_mdp_path:
        raise ValueError("PATH is required.")
    mdp_dir = os.path.join(base_path, "mdp")
    cmd = ["ansible-playbook",os.path.join(SHELL_DIR,"cp_dir.yml"),"-i","localhost,","-c","local","--extra-vars",json.dumps({"src":reference_mdp_path,"dst":mdp_dir})]
    try:
        subprocess.run(cmd)
        return True
    except:
        return False


@contextmanager
def _in_dir(p):
    o=os.getcwd();os.makedirs(p,exist_ok=True);os.chdir(p)
    try:yield p
    finally:os.chdir(o)

def fasta_to_pir(f,c):
    s=''.join(i.strip()for i in f.splitlines()if not i.startswith('>'))
    return f'>P1;{c}\nsequence:{c}:::::::0.00: 0.00\n{s}*\n'

def _best_dope_model(w,s):
    g=glob.glob(os.path.join(w,f'{s}.B9999*.pdb'))
    if not g:raise FileNotFoundError(f'モデル PDB が見つかりません: {os.path.join(w,f"{s}.B9999*.pdb")}')
    e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
    b,k=None,float('inf')
    for p in g:
        d=Selection(complete_pdb(e,p)).assess_dope()
        if d<k:k,b=d,p
    return b

def build_profile(workdir,query_pir,pdb_db_pir,query_code='TARGET',matrix_offset=-450,gap_penalties_1d=(-500,-50),max_aln_evalue=.01,n_prof_iterations=1):
    with _in_dir(workdir):
        log.verbose();e=Environ();s=SequenceDB(e);b='pdb_95.bin'
        if not os.path.exists(b):
            s.read(seq_database_file=pdb_db_pir,seq_database_format='PIR',chains_list='ALL',minmax_db_seq_len=(30,4000),clean_sequences=True)
            s.write(seq_database_file=b,seq_database_format='BINARY',chains_list='ALL')
        s.read(seq_database_file=b,seq_database_format='BINARY',chains_list='ALL')
        a=Alignment(e);a.append(file=query_pir,alignment_format='PIR',align_codes='ALL')
        p=a.to_profile()
        p.build(s,matrix_offset=matrix_offset,rr_file='${LIB}/blosum62.sim.mat',gap_penalties_1d=gap_penalties_1d,n_prof_iterations=n_prof_iterations,check_profile=False,max_aln_evalue=max_aln_evalue)
        p.write(file='build_profile.prf',profile_format='TEXT')
        a=p.to_alignment();a.write(file='build_profile.ali',alignment_format='PIR')
    return os.path.join(workdir,'build_profile.prf')

def compare_templates(workdir,template_list,atom_files_dir='.'):
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[atom_files_dir,'.'];a=Alignment(e)
        for p,c in template_list:
            a.append_model(Model(e,file=p,model_segment=(f'FIRST:{c}',f'LAST:{c}')),atom_files=p,align_codes=p+c)
        a.malign();a.malign3d();a.compare_structures();a.id_table(matrix_file='family.mat');e.dendrogram(matrix_file='family.mat',cluster_cut=-1.)
    return os.path.join(workdir,'compare.log')

def align2d_single(workdir,template_pdb,template_chain,query_pir,query_code,template_code=None,max_gap_length=50):
    b=os.path.splitext(os.path.basename(template_pdb))[0]
    if template_code is None:template_code=b+template_chain
    o=f'{query_code}-{template_code}.ali'
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb))or'.','.']
        a=Alignment(e)
        a.append_model(Model(e,file=template_pdb,model_segment=(f'FIRST:{template_chain}',f'LAST:{template_chain}')),align_codes=template_code,atom_files=template_pdb)
        a.append(file=query_pir,align_codes=query_code);a.align2d(max_gap_length=max_gap_length)
        a.write(file=o,alignment_format='PIR');a.write(file=o.replace('.ali','.pap'),alignment_format='PAP')
    return os.path.join(workdir,o)

def build_single_model(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def evaluate_model(workdir,pdb_file,output_profile=None,smoothing_window=15):
    b=os.path.splitext(os.path.basename(pdb_file))[0]
    if output_profile is None:output_profile=f'{b}.profile'
    with _in_dir(workdir):
        log.verbose();e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
        Selection(complete_pdb(e,os.path.abspath(pdb_file))).assess_dope(output='ENERGY_PROFILE NO_REPORT',file=output_profile,normalize_profile=True,smoothing_window=smoothing_window)
    return os.path.join(workdir,output_profile)

def salign_multiple_templates(workdir,template_list,atom_files_dir='.',output_ali='templates_mult.ali'):
    with _in_dir(workdir):
        log.verbose();e=Environ();e.io.atom_files_directory=[atom_files_dir,'.'];a=Alignment(e)
        for c,h in template_list:
            a.append_model(Model(e,file=c,model_segment=(f'FIRST:{h}',f'LAST:{h}')),atom_files=c,align_codes=c+h)
        for w,f,o in(((1.,0.,0.,0.,1.,0.),False,True),((1.,.5,1.,1.,1.,0.),False,True),((1.,1.,1.,1.,1.,0.),True,False)):
            a.salign(rms_cutoff=3.5,normalize_pp_scores=False,rr_file='$(LIB)/as1.sim.mat',overhang=30,gap_penalties_1d=(-450,-50),gap_penalties_3d=(0,3),gap_gap_score=0,gap_residue_score=0,dendrogram_file='templates.tree',alignment_type='tree',feature_weights=w,improve_alignment=True,fit=True,write_fit=f,write_whole_pdb=o,output='ALIGNMENT QUALITY')
        a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP');a.write(file=output_ali,alignment_format='PIR')
        a.salign(rms_cutoff=1.,normalize_pp_scores=False,rr_file='$(LIB)/as1.sim.mat',overhang=30,gap_penalties_1d=(-450,-50),gap_penalties_3d=(0,3),gap_gap_score=0,gap_residue_score=0,dendrogram_file='templates_quality.tree',alignment_type='progressive',feature_weights=[0]*6,improve_alignment=False,fit=False,write_fit=True,write_whole_pdb=False,output='QUALITY')
    return os.path.join(workdir,output_ali)

def align2d_multiple(workdir,templates_ali,query_pir,query_code,output_ali='query_mult.ali',max_gap_length=20):
    with _in_dir(workdir):
        log.verbose();e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib')
        a=Alignment(e);a.append(file=os.path.abspath(templates_ali),align_codes='all');b=len(a)
        a.append(file=os.path.abspath(query_pir),align_codes=query_code)
        a.salign(output='',max_gap_length=max_gap_length,gap_function=True,alignment_type='PAIRWISE',align_block=b,feature_weights=(1.,0.,0.,0.,0.,0.),overhang=0,gap_penalties_1d=(-450,0),gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),similarity_flag=True)
        a.write(file=output_ali,alignment_format='PIR');a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP')
    return os.path.join(workdir,output_ali)

def build_multi_model(workdir,ali_file,template_codes,sequence,n_models=5,assess_methods=None,use_hetatm=False):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.'];e.io.hetatm=use_hetatm
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def refine_loop(workdir,ini_model,sequence,loop_range,n_loop_models=10,md_level=None,atom_files_dir='.'):
    if md_level is None:md_level=refine.very_fast
    s,e=loop_range
    class _MyLoop(LoopModel):
        def select_loop_atoms(self):return Selection(self.residue_range(s,e))
    with _in_dir(workdir):
        log.verbose();v=Environ();v.io.atom_files_directory=[atom_files_dir,'.']
        m=_MyLoop(v,inimodel=os.path.abspath(ini_model),sequence=sequence)
        m.loop.starting_model=1;m.loop.ending_model=n_loop_models;m.loop.md_level=md_level;m.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.BL*.pdb')))

def build_model_with_ligand(workdir,ali_file,template_codes,sequence,restraint_atom_pairs=None,restraint_mean=3.5,restraint_stdev=.1,n_models=5):
    if restraint_atom_pairs is None:restraint_atom_pairs=[]
    a,m,s=restraint_atom_pairs,restraint_mean,restraint_stdev
    class _MyModel(AutoModel):
        def special_restraints(self,aln):
            r=self.restraints
            for i,j in a:r.add(forms.UpperBound(group=physical.upper_distance,feature=features.Distance(self.atoms[i],self.atoms[j]),mean=m,stdev=s))
    with _in_dir(workdir):
        e=Environ();e.io.hetatm=True;e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        x=_MyModel(e,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence)
        x.starting_model=1;x.ending_model=n_models;x.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def align2d_with_ss(workdir,template_pdb,template_chain,query_pir,query_code,output_ali=None,max_gap_length=50):
    b=os.path.splitext(os.path.basename(template_pdb))[0];t=b+template_chain
    if output_ali is None:output_ali=f'{query_code}-{t}.ali'
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb))or'.','.']
        a=Alignment(e)
        a.append_model(Model(e,file=template_pdb,model_segment=(f'FIRST:{template_chain}',f'LAST:{template_chain}')),align_codes=t,atom_files=template_pdb)
        a.append(file=os.path.abspath(query_pir),align_codes=query_code);a.align2d(max_gap_length=max_gap_length)
        a.write(file=output_ali,alignment_format='PIR')
        a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP',alignment_features='INDICES HELIX BETA')
    return os.path.join(workdir,output_ali)

def iterative_modeling(workdir,template_pdb,template_chain,query_pir,query_code,max_iterations=5,n_models=1):
    b=os.path.splitext(os.path.basename(template_pdb))[0];t=b+template_chain
    k,p=float('inf'),None
    for i in range(1,max_iterations+1):
        d=os.path.join(workdir,f'iter_{i:02d}');os.makedirs(d,exist_ok=True)
        print(f'\n=== 反復 {i}/{max_iterations} ===')
        a=align2d_with_ss(workdir=d,template_pdb=template_pdb,template_chain=template_chain,query_pir=query_pir,query_code=query_code)
        g=build_single_model(workdir=d,ali_file=a,template_code=t,sequence=query_code,n_models=n_models,assess_methods=(assess.DOPE,assess.GA341))
        if not g:
            print(f'  モデル生成に失敗しました (反復 {i})');break
        e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
        q,r=float('inf'),None
        for x in g:
            y=Selection(complete_pdb(e,x)).assess_dope()
            if y<q:q,r=y,x
        print(f'  最良 DOPE スコア: {q:.2f}  ({r})')
        if q<k:k,p=q,r
        else:
            print('  スコアが改善しませんでした。反復を終了します。');break
    print(f'\n最終ベストモデル: {p} (DOPE={k:.2f})')
    return p

def build_from_threading_ali(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    p=sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))
    e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
    s=[]
    for i in p:s.append((Selection(complete_pdb(e,i)).assess_dope(),i))
    s.sort()
    return[i for _,i in s]

def get_pdb(b,p):
    p=p.upper()
    w=os.path.abspath(os.path.join(b,"sys",p));os.makedirs(w,exist_ok=True)
    out_dir=os.path.abspath(os.path.join(b,"sys"))
    pdb=f"{w}/{p}.pdb";fa=f"{w}/{p}.fasta";pir=f"{w}/{p}.pir";ali=f"{w}/{p}.ali"
    clean=f"{out_dir}/clean.pdb"
    if not os.path.exists(pdb):r=requests.get(f"https://files.rcsb.org/download/{p}.pdb",timeout=60);r.raise_for_status();open(pdb,"wb").write(r.content)
    if not os.path.exists(fa):r=requests.get(f"https://www.rcsb.org/fasta/entry/{p}/download",timeout=60);r.raise_for_status();open(fa,"w").write(r.text)
    if not os.path.exists(pir):open(pir,"w").write(fasta_to_pir(open(fa).read(),"TARGET"))
    a=align2d_single(workdir=w,template_pdb=pdb,template_chain="A",query_pir=pir,query_code="TARGET")
    build_single_model(workdir=w,ali_file=a,template_code=f"{p}A",sequence="TARGET",n_models=1)
    b=_best_dope_model(w,"TARGET");shutil.copy2(b,clean);return clean
    
'''
def get_pdb(base_path, pdbid):
    if not base_path:
        raise ValueError("PATH is required.")
    if not pdbid:
        raise ValueError("pdbid is required.")
    sys_dir = os.path.join(base_path, "sys")
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    out_path = os.path.join(sys_dir, f"{pdbid}.pdb")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(r.content)
    return out_path
   
from shinka.llm.query import query
def simulation_set(base_path, simulation_information = None):    
    if not base_path:
        raise ValueError("PATH is required.")
    mdp_dir = os.path.join(base_path, "mdp")
    aimodel = "local/ai/gemma3-qat@http://172.17.0.1:12434/v1"
    cmds = [
        {"makemdp": "ions" , "prompts": "Generate a valid GROMACS ions.mdp file for ion addition (grompp before genion). Use typical settings for energy minimization preparation. Output only the mdp file content."},
        {"makemdp": "em" , "prompts": "Generate a valid GROMACS em.mdp file for energy minimization. Use steepest descent minimization and typical parameters used in protein simulations. Output only the mdp file content."},
        {"makemdp": "nvt" , "prompts": "Generate a GROMACS nvt.mdp file for NVT equilibration. Include temperature coupling (e.g., V-rescale thermostat) and position restraints for protein heavy atoms. Output only the mdp file content."},
        {"makemdp": "npt_br" , "prompts": "Generate a GROMACS npt.mdp file for NPT equilibration using the Berendsen pressure coupling method. This stage should stabilize the system density. Output only the mdp file content."},
        {"makemdp": "npt_pr" , "prompts": "Generate a GROMACS npt.mdp file using Parrinello-Rahman pressure coupling for proper NPT ensemble equilibration before production MD. Output only the mdp file content."},
        {"task_type": "watermodel", "prompts": "Guess the water model used in the simulation based on the provided information (e.g., TIP3P, SPC/E, TIP4P). Output only the water model name."},
        {"task_type": "force-field", "prompts": "Guess the force field used in the simulation based on the provided information (e.g., AMBER99SB-ILDN, CHARMM36, OPLS-AA). Output only the force field name."},
        {"task_type": "waterboxfile", "prompts": "Guess the pre-equilibrated water box filename used in the simulation (e.g., spc216.gro, tip3p.gro). Output only the filename."},
        {"task_type": "distance", "prompts": "Guess the minimum distance between the protein and the box edge in nm used in the simulation (e.g., 1.0, 1.5). Output only a number."}
    ]
    water_model, force_field,  waterboxfile, distance = None, None, None, None 
    try:
        for c in cmds:
            key = next(iter(c))
            if key == "makemdp":
                msg =f"Task:\n{c['prompts']}"
                if simulation_information:
                    msg = f"Simulation conditions from here:\n{simulation_information}\n\n{msg}"
                res = query(
                    model_name=aimodel,
                    msg=msg,
                    system_msg=("You are an expert in molecular dynamics simulations using GROMACS. Generate valid mdp configuration files. Follow the provided simulation conditions if available. Return only the mdp file content without explanations.")
                )
                content = getattr(res, "content", "")
                code = re.sub(r'^```.*?\n|```$', '', content, flags=re.MULTILINE).strip()
                with open(os.path.join(mdp_dir, f"{c['makemdp']}.mdp"), "w", encoding="utf-8") as f:
                    f.write(code)
            if key == "task_type":
                msg =f"Task:\n{c['prompts']}"
                if simulation_information:
                    msg = f"Simulation conditions from here:\n{simulation_information}\n{msg}"
                if distance:
                    msg = f"minimum distance between the protein and the box edge: {distance}\n{msg}"
                if waterboxfile:
                    msg = f"pre-equilibrated water box file: {waterboxfile}\n{msg}"
                if force_field:
                    msg = f"force field: {force_field}\n{msg}"
                if water_model:
                    msg = f"water model: {water_model}\n{msg}"
                res = query(
                    model_name=aimodel,
                    msg=msg,
                    system_msg="You are an information extraction assistant. Return only the requested value without explanations."
                )
                content = str(getattr(res, "content", "")).strip()
                if c["task_type"] == "watermodel":
                    water_model = content
                elif c["task_type"] == "force-field":
                    force_field = content
                elif c["task_type"] == "waterboxfile":
                    waterboxfile = content
                elif c["task_type"] == "distance":
                    try:
                        distance = float(content)
                    except ValueError:
                        distance = None
            print(f"Done {key}")
        return water_model, force_field,  waterboxfile, distance
    except:
        return False
'''


def system_build(base_path, pdb_path, FF, DISTANCE, WATER_MODEL, WATERBOXFILE, GMX="gmx", PYMOL="pymol"):
    if not base_path:
        raise ValueError("BASEPATH is required")
    if not pdb_path:
        raise ValueError("PDBPATH is required")
    if not FF:
        raise ValueError("forcefield is required")
    if not DISTANCE:
        raise ValueError("distance is required")
    if not WATER_MODEL:
        raise ValueError("water-model is required")
    if not WATERBOXFILE:
        raise ValueError("water box file is required")

    sys_dir     = os.path.join(base_path, "sys")
    mdp_dir     = os.path.join(base_path, "mdp")

    cle_pdb     = os.path.join(sys_dir, "clean.pdb")
    pro_gro     = os.path.join(sys_dir, "processed.gro")
    pro_top     = os.path.join(sys_dir, "processed.top")
    nbox_gro    = os.path.join(sys_dir, "newbox.gro")
    sol_top     = os.path.join(sys_dir, "solv.top")
    sol_gro     = os.path.join(sys_dir, "solv.gro")
    ions_tpr    = os.path.join(sys_dir, "ions.tpr")
    ions_mdp    = os.path.join(mdp_dir, "ions.mdp")
    ions_gro    = os.path.join(sys_dir, "ions.gro")
    ions_top    = os.path.join(sys_dir, "ions.top")
    md_gro      = os.path.join(sys_dir, "MD.gro")
    md_top      = os.path.join(sys_dir, "MD.top")

    cmds = [
        {"cmd": [PYMOL, "-cq", "-d", f"load {pdb_path}; remove not polymer; save {cle_pdb}; quit"]},
        {"cmd": [GMX, "pdb2gmx", "-f", cle_pdb, "-o", pro_gro, "-water", WATER_MODEL, "-p", pro_top, "-ff", FF]},
        {"cmd": [GMX, "editconf", "-f", pro_gro, "-o", nbox_gro, "-c", "-d", str(DISTANCE), "-bt", "cubic"]},
        {"cmd": ["ansible-playbook",os.path.join(SHELL_DIR,"cp_file.yml"),"-i","localhost,","-c","local","--extra-vars",json.dumps({"src":pro_top,"dst":sol_top})]},
        {"cmd": [GMX, "solvate", "-cp", nbox_gro, "-cs", WATERBOXFILE, "-o", sol_gro, "-p", sol_top]},
        {"cmd": [GMX, "grompp", "-f", ions_mdp, "-c", sol_gro, "-p", sol_top, "-o", ions_tpr, "-maxwarn", "1"]},
        {"cmd": ["ansible-playbook",os.path.join(SHELL_DIR,"cp_file.yml"),"-i","localhost,","-c","local","--extra-vars",json.dumps({"src":sol_top,"dst":ions_top})]},
        {"cmd-input": [GMX, "genion", "-s", ions_tpr, "-o", ions_gro, "-p", ions_top, "-pname", "NA", "-nname", "CL", "-neutral"], "input": "SOL\n"},
        {"cmd": ["ansible-playbook",os.path.join(SHELL_DIR,"cp_file.yml"),"-i","localhost,","-c","local","--extra-vars",json.dumps({"src":ions_gro,"dst":md_gro})]},
        {"cmd": ["ansible-playbook",os.path.join(SHELL_DIR,"cp_file.yml"),"-i","localhost,","-c","local","--extra-vars",json.dumps({"src":ions_top,"dst":md_top})]},
    ]
    try:
        for c in cmds:
            key = next(iter(c))
            if key == "cmd":
                subprocess.run(c["cmd"])
            elif key == "cmd-input":
                subprocess.run(c["cmd-input"], input=c["input"], text=True)
            else:
                return False
            print(f"Done {key}")
        return True
    except:
        return False


def minimization(base_path, GMX="gmx"):
    if not base_path:
        raise ValueError("PATH is required")
    sys_dir     = os.path.join(base_path, "sys")
    em_dir      = os.path.join(base_path, "em")
    mdp_dir     = os.path.join(base_path, "mdp")

    em_mdp      = os.path.join(mdp_dir, "em.mdp")
    md_gro      = os.path.join(sys_dir, "MD.gro")
    md_top      = os.path.join(sys_dir, "MD.top")
    em_tpr      = os.path.join(em_dir, "em.tpr")
    em_prefix   = os.path.join(em_dir, "em")

    cmds = [
        [GMX, "grompp", "-f", em_mdp, "-c", md_gro, "-p", md_top, "-o", em_tpr, "-maxwarn", "1"],
        [GMX, "mdrun","-v","-deffnm",em_prefix]
    ]
    try:
        for c in cmds:
            subprocess.run(c)
        return True
    except:
        return False

def nvt(base_path, GMX="gmx"):
    if not base_path:
        raise ValueError("PATH is required")
    sys_dir     = os.path.join(base_path, "sys")
    em_dir      = os.path.join(base_path, "em")
    nvt_dir     = os.path.join(base_path, "nvt")
    mdp_dir     = os.path.join(base_path, "mdp")

    md_top      = os.path.join(sys_dir, "MD.top")
    nvt_mdp     = os.path.join(mdp_dir, "nvt.mdp")
    em_gro      = os.path.join(em_dir, "em.gro")
    nvt_tpr     = os.path.join(nvt_dir, "nvt.tpr")
    nvt_prefix  = os.path.join(nvt_dir, "nvt")

    cmds = [
        [GMX, "grompp", "-f", nvt_mdp, "-c", em_gro, "-r", em_gro, "-p", md_top, "-o", nvt_tpr, "-maxwarn", "1"],
        [GMX, "mdrun", "-deffnm", nvt_prefix]
    ]
    try:
        for c in cmds:
            subprocess.run(c)
        return True
    except:
        return False

def npt_br(base_path, GMX="gmx"):
    if not base_path:
        raise ValueError("PATH is required")
    sys_dir     = os.path.join(base_path, "sys")
    nvt_dir     = os.path.join(base_path, "nvt")
    npt_br_dir  = os.path.join(base_path, "npt_br")
    mdp_dir     = os.path.join(base_path, "mdp")
    md_top      = os.path.join(sys_dir, "MD.top")
    npt_br_mdp  = os.path.join(mdp_dir, "npt_br.mdp")
    nvt_gro     = os.path.join(nvt_dir, "nvt.gro")
    npt_br_tpr  = os.path.join(npt_br_dir, "npt_br.tpr")
    npt_br_prefix = os.path.join(npt_br_dir, "npt_br")
    cmds = [
        [GMX, "grompp", "-f", npt_br_mdp, "-c", nvt_gro, "-r", nvt_gro, "-p", md_top, "-o", npt_br_tpr, "-maxwarn", "1"],
        [GMX, "mdrun", "-deffnm", npt_br_prefix]
    ]
    try:
        for c in cmds:
            subprocess.run(c, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except Exception:
        return False
    

def npt_pr(base_path, GMX="gmx"):
    if not base_path:
        raise ValueError("PATH is required")

    sys_dir     = os.path.join(base_path, "sys")
    npt_br_dir  = os.path.join(base_path, "npt_br")
    npt_pr_dir  = os.path.join(base_path, "npt_pr")
    mdp_dir     = os.path.join(base_path, "mdp")

    md_top      = os.path.join(sys_dir, "MD.top")
    npt_pr_mdp  = os.path.join(mdp_dir, "npt_pr.mdp")
    npt_br_gro  = os.path.join(npt_br_dir, "npt_br.gro")

    npt_pr_tpr  = os.path.join(npt_pr_dir, "npt_pr.tpr")
    npt_pr_prefix = os.path.join(npt_pr_dir, "npt_pr")

    cmds = [
        [GMX, "grompp", "-f", npt_pr_mdp, "-c", npt_br_gro, "-r", npt_br_gro, "-p", md_top, "-o", npt_pr_tpr, "-maxwarn", "1"],
        [GMX, "mdrun", "-deffnm", npt_pr_prefix]
    ]

    try:
        for c in cmds:
            subprocess.run(c, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except Exception:
        return False

def md(base_path, GMX="gmx"):
    if not base_path:
        raise ValueError("PATH is required")

    sys_dir     = os.path.join(base_path, "sys")
    npt_pr_dir  = os.path.join(base_path, "npt_pr")
    md_dir      = os.path.join(base_path, "md")
    mdp_dir     = os.path.join(base_path, "mdp")

    md_top      = os.path.join(sys_dir, "MD.top")
    md_mdp      = os.path.join(mdp_dir, "md.mdp")
    npt_pr_gro  = os.path.join(npt_pr_dir, "npt_pr.gro")

    md_tpr      = os.path.join(md_dir, "md.tpr")
    md_prefix   = os.path.join(md_dir, "md")

    cmds = [
        [GMX, "grompp", "-f", md_mdp, "-c", npt_pr_gro, "-p", md_top, "-o", md_tpr],
        [GMX, "mdrun", "-deffnm", md_prefix]
    ]

    try:
        for c in cmds:
            subprocess.run(c, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    except Exception:
        return False
