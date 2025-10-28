from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple
import os, random

from .base import Engine
from ..gnina_wrapper import dock as gnina_dock, rescore as gnina_rescore
from ..pdbqt_io import to_pdbqt_receptor, to_pdbqt_ligand

WORK_DIR = Path("data/work/docking")

class GninaEngine(Engine):
    name = "gnina"

    def dock(self, protein_pdb_in: str|Path, ligand_sdf_in: str|Path, cfg: Dict[str,Any]) -> List[str]:
        WORK_DIR.mkdir(parents=True, exist_ok=True)
        protein_pdb_in = Path(protein_pdb_in)
        ligand_sdf_in  = Path(ligand_sdf_in)

        # Prep receptor/ligand PDBQT
        receptor_pdbqt = WORK_DIR / f"{protein_pdb_in.stem}.pdbqt"
        if not receptor_pdbqt.exists():
            to_pdbqt_receptor(protein_pdb_in, receptor_pdbqt)
        ligand_pdbqt = WORK_DIR / f"{ligand_sdf_in.stem}.pdbqt"
        if not ligand_pdbqt.exists():
            to_pdbqt_ligand(ligand_sdf_in, ligand_pdbqt)

        dcfg = (cfg.get("docking") or {})
        gcfg = (dcfg.get("gnina") or {})
        box  = (dcfg.get("box") or {})
        center = tuple(box.get("center") or (0.0,0.0,0.0))  # Expect tests to have validated box center
        size   = box.get("size") or (20.0,20.0,20.0)
        if isinstance(size, (int,float)): size = (float(size),)*3
        num_modes = int(dcfg.get("poses", 20))
        exhaust   = int(dcfg.get("exhaustiveness", 8))
        base_seed = int(dcfg.get("seed", 42))
        num_seeds = int(dcfg.get("num_seeds", 1))
        threads   = int((cfg.get("runtime") or {}).get("cpu_threads", max(1, (os.cpu_count() or 2)//2)))
        device    = str(gcfg.get("device", "auto"))
        model     = gcfg.get("model") or "default"

        out_paths: List[str] = []
        for i in range(num_seeds):
            seed = base_seed if num_seeds == 1 else random.randint(1, 10_000_000)
            out_sdf = WORK_DIR / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}_g{seed}_poses.sdf"
            log     = WORK_DIR / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}_g{seed}_gnina.log"
            res = gnina_dock(
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=ligand_pdbqt,
                center=(float(center[0]), float(center[1]), float(center[2])),
                size=(float(size[0]), float(size[1]), float(size[2])),
                num_modes=num_modes,
                exhaustiveness=exhaust,
                seed=seed,
                cpu_threads=threads,
                device=device,
                model=model,
                out_sdf=out_sdf,
                log_path=log,
            )
            out_paths.append(out_sdf.as_posix())
        return out_paths

    def rescore(self, poses_sdf: str|Path, receptor_pdbqt: str|Path, cfg: Dict[str,Any]) -> list[dict]:
        dcfg = (cfg.get("docking") or {})
        gcfg = (dcfg.get("gnina") or {})
        device = str(gcfg.get("device", "auto"))
        model  = gcfg.get("model") or "default"
        out_sdf = Path(poses_sdf).with_name(Path(poses_sdf).stem + "_gnina_rescored.sdf")
        log     = out_sdf.with_suffix(".log")
        return gnina_rescore(Path(poses_sdf), Path(receptor_pdbqt), device=device, model=model, out_sdf=out_sdf, log_path=log)
