"""
UNIFIED ENVIRONMENT MANAGER
Replaces scattered shell scripts with single Python tool
"""
import argparse
import subprocess
import sys
from pathlib import Path

ENV_NAME = "ovmpk"
CONDA_CMD = "mamba"  # or "conda"

def create_env():
    """Create new conda environment"""
    cmd = [
        CONDA_CMD, "env", "create",
        "-n", ENV_NAME,
        "--file", "environment.yaml"
    ]
    subprocess.run(cmd, check=True)

def install_deps():
    """Install system dependencies"""
    deps = ["cuda-toolkit", "openmm"]
    cmd = ["sudo", "apt-get", "install", "-y"] + deps
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    subparsers.add_parser('create', help='Create new environment')
    subparsers.add_parser('install', help='Install system deps')
    subparsers.add_parser('freeze', help='Freeze current environment')
    
    args = parser.parse_args()
    
    if args.command == "create":
        create_env()
    elif args.command == "install":
        install_deps()
    elif args.command == "freeze":
        subprocess.run([CONDA_CMD, "env", "export", "-n", ENV_NAME], check=True)

if __name__ == "__main__":
    main()
