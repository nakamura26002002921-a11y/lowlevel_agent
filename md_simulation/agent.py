import subprocess
import os
#from shinka.llm.query import query
import re
import requests
import shlex
import shutil
from pathlib import Path
import json

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

'''
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
        [GMX,"grompp","-f",em_mdp,"-c",md_gro,"-p",md_top,"-o",em_tpr],
        [GMX,"mdrun","-v","-deffnm",em_prefix]
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
        [GMX, "grompp", "-f", nvt_mdp, "-c", em_gro, "-r", em_gro, "-p", md_top, "-o", nvt_tpr],
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
        [GMX, "grompp", "-f", npt_br_mdp, "-c", nvt_gro, "-r", nvt_gro, "-p", md_top, "-o", npt_br_tpr],
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
