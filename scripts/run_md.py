"""
UNIFIED MD WORKFLOW
Combines md_gpu_prepare.py and md_gpu_run.py
"""
import argparse
from openmm import app
import logging

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Preparation subcommand
    prep_parser = subparsers.add_parser('prepare')
    prep_parser.add_argument("input", help="Input structure file")
    prep_parser.add_argument("-o", "--output", default="system.xml")
    
    # Run subcommand  
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument("system", help="System XML file")
    run_parser.add_argument("-t", "--time", type=float, default=10.0)
    
    args = parser.parse_args()
    
    if args.command == "prepare":
        from md_prep import prepare_system
        prepare_system(args.input, args.output)
    elif args.command == "run":
        from md_run import simulate_system
        simulate_system(args.system, args.time)

if __name__ == "__main__":
    main()
