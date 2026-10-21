# Virtual-screening-suite (VS-Suite)
VS-Suite, an In-house, integrated virtual screening pipeline that orchestrates three consensus docking approaches (1) Boltz-2, co-folding model (2) AutoDock Vina, empirical scoring function; and (3) GNINA, which augments Vina-style docking with  CNN rescoring within a automated workflow accessible through a browser-based GUI. VS-Suite integrates three independent molecular docking methods into one pipeline with a browser-based GUI:

| Tool | Role |
|------|------|
| **Boltz-2** | AI-based protein-ligand co-folding + binding affinity (IC₅₀, P(binder)) |
| **AutoDock Vina** | Classical empirical docking (ΔG kcal/mol) |
| **GNINA** | Vina + deep-learned CNN scoring (pose score + pK) |

---

##  Key Features

- **Input** - protein FASTA + ligands CSV → complete docking results
- **Browser GUI** - upload files, click Run, view results without command line
- **Cross-tool consensus ranking** - ranks hits by agreement across all 3 tools
- **Bulk downloads** - top-10 (all poses) and top-20 (best pose per tool) as ZIP
- **MD-ready outputs** - SDF poses ready for GROMACS/AMBER parameterization
- **2D interaction analysis** - protein-ligand interactions via PLIP
- **Validated** - tested on 58 antiviral compounds against SARS-CoV-2 Mpro

---

## Let's Start

```bash
# 1. Clone
git clone https://github.com/Vsubuchem-lab/Virtual-screening-suite-VS-Suite.git
cd Virtual-screening-suite-VS-Suite-main

# 2. Install (one-time ~30 min)
bash install.sh

# 3. Launch GUI
conda activate autodock-main
streamlit run gui/VS_GUI.py
# Open http://localhost:8501
```

---

## Installation — Full Guide

### Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| OS | Linux Ubuntu 20.04+ | Ubuntu 22.04 |
| Python | 3.10 | 3.10 |
| RAM | 16 GB | 32 GB |
| GPU | None (CPU fallback) | NVIDIA 12 GB+ VRAM |
| CUDA | — | 11.8 or 12.x |
| Disk | 50 GB | 200 GB |

> **Windows users**: use WSL2 with Ubuntu.  
> **macOS**: not supported (no CUDA for Boltz-2/GNINA).

---

### Step 1 — Install Conda

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# Close and reopen terminal
conda --version    # confirm
```

---

### Step 2 — Clone VS-Suite

```bash
git clone https://github.com/Vsubuchem-lab/Virtual-screening-suite-VS-Suite.git
cd VS-Suite-GUI-Wrapper
```

---

### Step 3 — Run the Installer

```bash
bash install.sh
```

This creates three conda environments:

| Env | Contents | Time |
|-----|----------|------|
| `autodock-main` | Python 3.10, RDKit, Streamlit, Biopython | ~3-5 min |
| `autodock-prep` | OpenBabel, PDB2PQR, PROPKA, Meeko, PLIP | ~10-15 min |
| `boltz2` | Boltz-2 + PyTorch | ~5 min + 5 GB model download |

---

### Step 4 — Install AutoDock Vina (manual binary)

```bash
mkdir -p ~/tools && cd ~/tools
wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
tar -xzf autodock_vina_1_1_2_linux_x86.tgz
echo 'export PATH="$HOME/tools/autodock_vina_1_1_2_linux_x86/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
vina --version    # confirm
```

---

### Step 5 — Install GNINA (manual binary)

```bash
mkdir -p ~/tools/gnina && cd ~/tools/gnina
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
chmod +x gnina
echo 'export PATH="$HOME/tools/gnina:$PATH"' >> ~/.bashrc
source ~/.bashrc
gnina --version    # confirm
```

---

### Step 6 — Verify Everything Works

```bash
cd ~/path/to/VS-Suite-GUI-Wrapper
conda activate autodock-main
python tests/test_phase1.py
```
