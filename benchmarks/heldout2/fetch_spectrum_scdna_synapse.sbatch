#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_spectrum_scdna
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Claim 4 primary truth: MSK SPECTRUM DLP+ scDNA (Synapse syn66366960). Requires ~/.synapseConfig (see the .py
# docstring). Not submitted until the user has set up the token. BLIND: download and checksum only.
# python without -I: the synapse env takes idna from ~/.local (user site), which -I hides. cwd is the repo, not the
# download dir, so no downloaded file sits on the import path.
set -euo pipefail
micromamba run -r /gpfs/commons/home/jameslee/micromamba -n synapse \
  python -B benchmarks/heldout2/fetch_spectrum_scdna_synapse.py data/external/spectrum_ov_scdna_synapse
