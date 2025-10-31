"""
UNIFIED PIPELINE RUNNER v2
Combines all pipeline variants with backward compatibility
"""
import argparse
import logging
import sys
from pathlib import Path

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Configuration file")
    parser.add_argument("--simple", action="store_true", 
                      help="Simplified execution mode")
    parser.add_argument("--log", default="INFO", 
                      choices=["DEBUG","INFO","WARNING","ERROR"],
                      help="Logging level")
    args = parser.parse_args()
    
    logging.basicConfig(level=args.log,
                      format='%(asctime)s [%(levelname)s] %(message)s')
    
    # Core functionality remains unchanged
    from pipeline.core import Pipeline
    pipeline = Pipeline(args.config)
    
    if args.simple:
        pipeline.run_simple()
    else:
        pipeline.run_full()

if __name__ == "__main__":
    run()
