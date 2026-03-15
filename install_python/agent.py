import subprocess
from pathlib import Path
BASE_DIR = Path(__file__).parent
def installpython(version, path):
    playbook_path = BASE_DIR / "install_python.yml"
    cmd = [
        "ansible-playbook",
        "-i", "hosts",
        str(playbook_path),
        "-e", f"VERSION={version} INSTALL_PATH={path}"
    ]
    print("[python_agent] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)