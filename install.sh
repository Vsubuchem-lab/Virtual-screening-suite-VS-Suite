#!/usr/bin/env bash
# VS-Suite installer
# Usage: bash install.sh

set -e

echo "========================================================"
echo "  VS-Suite Installer"
echo "  Centre for Advanced Computation and Theoretical"
echo "  Sciences (ACTS), SRMIST"
echo "========================================================"
echo ""

# Check conda
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found."
    echo "Install Miniforge: https://github.com/conda-forge/miniforge"
    exit 1
fi
echo "conda: $(conda --version)"

# Speed up solver
conda install -n base conda-libmamba-solver -y -q 2>/dev/null || true
conda config --set solver libmamba 2>/dev/null || true

# Env 1: orchestrator
if conda env list | grep -q "^autodock-main "; then
    echo "✓ autodock-main already exists"
else
    echo "→ Creating autodock-main..."
    conda env create -f envs/main.yml
fi

# Env 2: prep tools
if conda env list | grep -q "^autodock-prep "; then
    echo "✓ autodock-prep already exists"
else
    echo "→ Creating autodock-prep..."
    conda env create -f envs/prep.yml
    conda run -n autodock-prep pip install scipy gemmi prody pydantic
    conda run -n autodock-prep conda install -c conda-forge plip -y 2>/dev/null || \
        conda run -n autodock-prep pip install plip 2>/dev/null || true
fi

# Env 3: Boltz-2
if conda env list | grep -q "^boltz2 "; then
    echo "✓ boltz2 already exists"
else
    echo "→ Creating boltz2..."
    conda create -n boltz2 python=3.10 -y -q
    conda run -n boltz2 pip install boltz torch
fi

# Install VS-Suite package
echo "→ Installing VS-Suite..."
conda run -n autodock-main pip install -e . -q

# Check binaries
echo ""
echo "→ Checking Vina + GNINA..."
VINA_OK=0; GNINA_OK=0
command -v vina  >/dev/null 2>&1 && { echo "✓ vina:  $(which vina)";  VINA_OK=1;  } || echo "✗ vina  NOT found"
command -v gnina >/dev/null 2>&1 && { echo "✓ gnina: $(which gnina)"; GNINA_OK=1; } || echo "✗ gnina NOT found"

echo ""
echo "========================================================"
echo "  Conda environments ready ✓"
echo "========================================================"

if [ $VINA_OK -eq 0 ] || [ $GNINA_OK -eq 0 ]; then
    echo ""
    echo "ACTION REQUIRED — install missing binaries:"
    [ $VINA_OK -eq 0 ] && cat << 'VINA'

  AutoDock Vina:
    mkdir -p ~/tools && cd ~/tools
    wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
    tar -xzf autodock_vina_1_1_2_linux_x86.tgz
    echo 'export PATH="$HOME/tools/autodock_vina_1_1_2_linux_x86/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
VINA
    [ $GNINA_OK -eq 0 ] && cat << 'GNINA'

  GNINA:
    mkdir -p ~/tools/gnina && cd ~/tools/gnina
    wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
    chmod +x gnina
    echo 'export PATH="$HOME/tools/gnina:$PATH"' >> ~/.bashrc
    source ~/.bashrc
GNINA
fi

echo ""
echo "  Next steps:"
echo "    conda activate autodock-main"
echo "    python tests/test_phase1.py    # smoke test"
echo "    streamlit run gui/VS_GUI.py    # launch GUI"
echo "    # open http://localhost:8501"
echo "========================================================"
