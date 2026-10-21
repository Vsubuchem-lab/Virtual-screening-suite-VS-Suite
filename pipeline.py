"""End-to-end docking pipeline orchestrator.

Single function: protein FASTA + ligands CSV -> all 3 tools, all ligands.

Usage:
    from core.pipeline import run_pipeline
    result = run_pipeline(
        fasta="examples/test_protein.fasta",
        ligands_csv="examples/test_ligands.csv",
        output_dir="results/run1",
        chain_mode="monomer",
        use_msa=True,
        run_vina=True,
        run_gnina=True,
    )
"""
import json
import time
import sys
import traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from core.utils.logger import setup_logger
from core.boltz.yaml_builder import build_boltz_yamls
from core.boltz.runner import run_boltz
from core.boltz.parser import parse_boltz_output
from core.prep.protonate import protonate_protein
from core.prep.protein_pdbqt import protein_to_pdbqt
from core.prep.ligand_pdbqt import ligand_to_pdbqt
from core.docking.vina import run_vina
from core.docking.gnina import run_gnina


def _step(label, fn, *args, **kwargs):
    """Wrap a step with timing + error capture."""
    start = time.time()
    try:
        result = fn(*args, **kwargs)
        return {"ok": True, "elapsed": round(time.time() - start, 1), "result": result}
    except Exception as e:
        return {
            "ok": False,
            "elapsed": round(time.time() - start, 1),
            "error": str(e),
            "trace": traceback.format_exc(),
        }


def run_one_ligand(yaml_info,
                   project_dir,
                   ligand_chain_id,
                   use_msa: bool = True,
                   use_potentials: bool = True,
                   diffusion_samples: int = 1,
                   recycling_steps: int = 3,
                   no_kernels: bool = False,
                   run_vina_flag: bool = True,
                   run_gnina_flag: bool = True,
                   vina_exhaustiveness: int = 16,
                   gnina_exhaustiveness: int = 16,
                   gnina_use_gpu: bool = True,
                   log=None) -> dict:
    """Process one ligand end-to-end."""
    title = yaml_info["title"]
    safe_title = yaml_info["safe_title"]
    yaml_path = yaml_info["yaml_path"]

    if log:
        log.info(f"\n{'='*70}\n  LIGAND: {title}\n{'='*70}")

    project_dir = Path(project_dir)
    boltz_out = project_dir / "boltz_runs"
    parsed_dir = project_dir / "boltz_parsed" / safe_title
    docking_dir = project_dir / "docking" / safe_title

    summary = {
        "title": title,
        "safe_title": safe_title,
        "smiles": yaml_info.get("smiles"),
        "n_atoms_heavy": yaml_info.get("n_atoms_heavy"),
        "yaml": yaml_path,
        "stages": {},
        "results": {
            "boltz": None,
            "vina": None,
            "gnina": None,
        },
    }

    # ----- STAGE 1: Boltz-2 -----
    if log:
        log.info("[Stage 1/4] Running Boltz-2...")
    s = _step("boltz_run",
              run_boltz,
              yaml_path=yaml_path,
              output_dir=boltz_out,
              use_msa_server=use_msa,
              use_potentials=use_potentials,
              output_format="pdb",
              diffusion_samples=diffusion_samples,
              recycling_steps=recycling_steps,
              no_kernels=no_kernels,
              stream=False,
              timeout=3600)
    summary["stages"]["boltz_run"] = {k: v for k, v in s.items() if k != "result"}
    if not s["ok"] or not s["result"].get("success"):
        if log:
            log.error(f"  Boltz failed: {s.get('error') or s['result'].get('error')}")
        return summary
    boltz_run_info = s["result"]
    if log:
        log.info(f"  Boltz done in {s['elapsed']}s")

    # ----- STAGE 2: Parse Boltz output -----
    if log:
        log.info("[Stage 2/4] Parsing Boltz output...")
    s = _step("boltz_parse",
              parse_boltz_output,
              predictions_dir=boltz_run_info["predictions_dir"],
              output_dir=parsed_dir,
              ligand_chain_id=ligand_chain_id)
    summary["stages"]["boltz_parse"] = {k: v for k, v in s.items() if k != "result"}
    if not s["ok"]:
        if log:
            log.error(f"  Parse failed: {s.get('error')}")
        return summary
    parsed = s["result"]
    summary["results"]["boltz"] = {
        "complex_pdb": parsed["complex_pdb"],
        "protein_pdb": parsed["protein_pdb"],
        "ligand_pdb": parsed["ligand_pdb"],
        "ligand_sdf": parsed["ligand_sdf"],
        "box": parsed["box"],
        "confidence": parsed["confidence"],
        "affinity": parsed["affinity"],
        "summary_path": parsed["summary_path"],
    }
    if log:
        aff = parsed.get("affinity", {})
        conf = parsed.get("confidence", {})
        log.info(f"  ipTM={conf.get('iptm', '-')}  pLDDT={conf.get('complex_plddt', '-')}")
        log.info(f"  Boltz IC50={aff.get('ic50_uM', '-')}uM  P_binder={aff.get('affinity_probability_binary', '-')}")

    # ----- STAGE 3: Prep for docking (only if doing Vina or GNINA) -----
    if not (run_vina_flag or run_gnina_flag):
        return summary

    if log:
        log.info("[Stage 3/4] Preparing PDBQT files for docking...")
    docking_dir.mkdir(parents=True, exist_ok=True)

    # Protein: protonate + PDBQT
    s = _step("protonate",
              protonate_protein,
              input_pdb=parsed["protein_pdb"],
              output_dir=docking_dir,
              ph=7.4,
              use_propka=True)
    summary["stages"]["protonate"] = {k: v for k, v in s.items() if k != "result"}
    if not s["ok"]:
        if log:
            log.error(f"  Protonate failed: {s.get('error')}")
        return summary
    protonated_pdb = s["result"]["protonated_pdb"]

    receptor_pdbqt = docking_dir / "receptor.pdbqt"
    s = _step("receptor_pdbqt",
              protein_to_pdbqt,
              input_pdb=protonated_pdb,
              output_pdbqt=receptor_pdbqt,
              charge_method="gasteiger")
    summary["stages"]["receptor_pdbqt"] = {k: v for k, v in s.items() if k != "result"}
    if not s["ok"]:
        if log:
            log.error(f"  Receptor PDBQT failed: {s.get('error')}")
        return summary

    # Ligand SDF -> PDBQT
    ligand_pdbqt = docking_dir / "ligand.pdbqt"
    s = _step("ligand_pdbqt",
              ligand_to_pdbqt,
              input_path=parsed["ligand_sdf"],
              output_pdbqt=ligand_pdbqt)
    summary["stages"]["ligand_pdbqt"] = {k: v for k, v in s.items() if k != "result"}
    if not s["ok"]:
        if log:
            log.error(f"  Ligand PDBQT failed: {s.get('error')}")
        return summary
    if log:
        log.info("  PDBQTs ready")

    box = parsed["box"]
    center = box["center"]
    size = box["size"]

    # ----- STAGE 4a: Vina -----
    if run_vina_flag:
        if log:
            log.info("[Stage 4a/4] Running Vina...")
        s = _step("vina",
                  run_vina,
                  receptor_pdbqt=str(receptor_pdbqt),
                  ligand_pdbqt=str(ligand_pdbqt),
                  output_dir=docking_dir / "vina",
                  center=center,
                  size=size,
                  exhaustiveness=vina_exhaustiveness,
                  n_poses=9,
                  cpu=4,
                  tag=f"{safe_title}_vina")
        summary["stages"]["vina"] = {k: v for k, v in s.items() if k != "result"}
        if s["ok"] and s["result"].get("success"):
            r = s["result"]
            summary["results"]["vina"] = {
                "best_affinity_kcal_mol": r["best_affinity_kcal_mol"],
                "n_poses": r["n_poses"],
                "output_pdbqt": r["output_pdbqt"],
                "poses": r["poses"],
            }
            if log:
                log.info(f"  Vina best ΔG = {r['best_affinity_kcal_mol']} kcal/mol")
        else:
            if log:
                log.error(f"  Vina failed: {s.get('error') or s.get('result', {}).get('error')}")

    # ----- STAGE 4b: GNINA -----
    if run_gnina_flag:
        if log:
            log.info("[Stage 4b/4] Running GNINA...")
        s = _step("gnina",
                  run_gnina,
                  receptor_pdbqt=str(receptor_pdbqt),
                  ligand_pdbqt=str(ligand_pdbqt),
                  output_dir=docking_dir / "gnina",
                  center=center,
                  size=size,
                  cnn_scoring="rescore",
                  exhaustiveness=gnina_exhaustiveness,
                  n_poses=9,
                  cpu=4,
                  use_gpu=gnina_use_gpu,
                  tag=f"{safe_title}_gnina")
        summary["stages"]["gnina"] = {k: v for k, v in s.items() if k != "result"}
        if s["ok"] and s["result"].get("success"):
            r = s["result"]
            summary["results"]["gnina"] = {
                "best_affinity_kcal_mol": r["best_affinity_kcal_mol"],
                "best_cnn_pose_score": r.get("best_cnn_pose_score"),
                "best_cnn_affinity_pK": r.get("best_cnn_affinity_pK"),
                "n_poses": r["n_poses"],
                "output_sdf": r["output_sdf"],
                "poses": r["poses"],
            }
            if log:
                log.info(f"  GNINA best ΔG = {r['best_affinity_kcal_mol']} kcal/mol, CNN_pK = {r.get('best_cnn_affinity_pK')}")
        else:
            if log:
                log.error(f"  GNINA failed: {s.get('error') or s.get('result', {}).get('error')}")

    return summary


def run_pipeline(fasta,
                 ligands_csv,
                 output_dir,
                 chain_mode: str = "monomer",
                 use_msa: bool = True,
                 use_potentials: bool = True,
                 diffusion_samples: int = 1,
                 recycling_steps: int = 3,
                 no_kernels: bool = False,
                 run_vina: bool = True,
                 run_gnina: bool = True,
                 vina_exhaustiveness: int = 16,
                 gnina_exhaustiveness: int = 16,
                 gnina_use_gpu: bool = True,
                 compute_affinity: bool = True,
                 affinity_max_atoms: int = 128) -> dict:
    """End-to-end docking pipeline.

    Args:
        fasta: protein FASTA file
        ligands_csv: CSV (title, smiles)
        output_dir: where everything goes
        chain_mode: monomer | dimer | all
        use_msa: use MSA server in Boltz
        use_potentials: Boltz inference-time potentials
        run_vina, run_gnina: toggle docking tools
        compute_affinity: run Boltz affinity prediction
        affinity_max_atoms: skip affinity if ligand bigger

    Returns:
        master dict with per-ligand results and overall summary
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logger(log_dir=output_dir / "logs", log_name="pipeline")
    log.info(f"\n{'#'*70}\n  AUTO-DOCK PIPELINE\n{'#'*70}")
    log.info(f"FASTA:    {fasta}")
    log.info(f"Ligands:  {ligands_csv}")
    log.info(f"Output:   {output_dir}")
    log.info(f"Mode:     {chain_mode}, MSA={use_msa}, Vina={run_vina}, GNINA={run_gnina}")

    # Step A: Build Boltz YAMLs
    log.info("\n[A] Building Boltz YAMLs...")
    boltz_inputs_dir = output_dir / "boltz_inputs"
    yaml_infos = build_boltz_yamls(
        fasta_path=fasta,
        ligands_csv=ligands_csv,
        output_dir=boltz_inputs_dir,
        chain_mode=chain_mode,
        compute_affinity=compute_affinity,
        affinity_max_atoms=affinity_max_atoms,
    )
    log.info(f"  Built {len(yaml_infos)} YAMLs")

    # Step B: Loop ligands
    pipeline_start = time.time()
    per_ligand = []
    for i, yi in enumerate(yaml_infos, 1):
        log.info(f"\n[B {i}/{len(yaml_infos)}] Processing: {yi['title']}")
        result = run_one_ligand(
            yaml_info=yi,
            project_dir=output_dir,
            ligand_chain_id=yi["ligand_chain_id"],
            use_msa=use_msa,
            use_potentials=use_potentials,
            diffusion_samples=diffusion_samples,
            recycling_steps=recycling_steps,
            no_kernels=no_kernels,
            run_vina_flag=run_vina,
            run_gnina_flag=run_gnina,
            vina_exhaustiveness=vina_exhaustiveness,
            gnina_exhaustiveness=gnina_exhaustiveness,
            gnina_use_gpu=gnina_use_gpu,
            log=log,
        )
        per_ligand.append(result)
        # Save intermediate after each ligand
        master = output_dir / "results_master.json"
        master.write_text(json.dumps({
            "fasta": str(fasta),
            "ligands_csv": str(ligands_csv),
            "chain_mode": chain_mode,
            "n_ligands": len(yaml_infos),
            "n_completed": i,
            "elapsed_seconds": round(time.time() - pipeline_start, 1),
            "per_ligand": per_ligand,
        }, indent=2, default=str))

    elapsed = time.time() - pipeline_start
    log.info(f"\n{'#'*70}\n  PIPELINE DONE in {elapsed/60:.1f} min\n{'#'*70}")

    # Compact comparison table
    log.info("\nCOMPARISON TABLE")
    log.info(f"{'Ligand':<20} {'Boltz IC50uM':>15} {'P_bind':>8} {'Vina dG':>10} {'GNINA dG':>10} {'GNINA CNNpK':>12}")
    log.info("-" * 80)
    for L in per_ligand:
        b = L["results"].get("boltz") or {}
        v = L["results"].get("vina") or {}
        g = L["results"].get("gnina") or {}
        aff = b.get("affinity", {})
        log.info(
            f"{L['title']:<20} "
            f"{str(aff.get('ic50_uM', '-')):>15} "
            f"{str(aff.get('affinity_probability_binary', '-'))[:6]:>8} "
            f"{str(v.get('best_affinity_kcal_mol', '-')):>10} "
            f"{str(g.get('best_affinity_kcal_mol', '-')):>10} "
            f"{str(g.get('best_cnn_affinity_pK', '-')):>12}"
        )

    return {
        "fasta": str(fasta),
        "ligands_csv": str(ligands_csv),
        "output_dir": str(output_dir),
        "chain_mode": chain_mode,
        "n_ligands": len(yaml_infos),
        "elapsed_seconds": round(elapsed, 1),
        "per_ligand": per_ligand,
        "master_json": str(output_dir / "results_master.json"),
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--ligands", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--chain-mode", default="monomer", choices=["monomer", "dimer", "all"])
    parser.add_argument("--no-msa", action="store_true")
    parser.add_argument("--no-vina", action="store_true")
    parser.add_argument("--no-gnina", action="store_true")
    parser.add_argument("--no-kernels", action="store_true")
    args = parser.parse_args()

    result = run_pipeline(
        fasta=args.fasta,
        ligands_csv=args.ligands,
        output_dir=args.out,
        chain_mode=args.chain_mode,
        use_msa=not args.no_msa,
        run_vina=not args.no_vina,
        run_gnina=not args.no_gnina,
        no_kernels=args.no_kernels,
    )
    print(f"\nMaster JSON: {result['master_json']}")
