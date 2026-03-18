#!/usr/bin/env bash
set -euo pipefail

# Install three tools into separate conda environments.
# Rationale: MobileElementFinder is distributed via pip, whereas
# IntegronFinder and ISEScan are available via Bioconda.
# Separate envs avoid dependency conflicts and make updates easier.

ME_ENV="mefinder_env"
IF_ENV="integronfinder_env"
IS_ENV="isescan_env"

# Configure channels (safe to run multiple times)
conda config --add channels bioconda || true
conda config --add channels conda-forge || true
conda config --set channel_priority strict || true

# 1) MobileElementFinder env
# BLAST is required; KMA is optional and omitted here.
conda create -y -n "${ME_ENV}" python=3.10 blast pip
conda run -n "${ME_ENV}" python -m pip install --upgrade pip
conda run -n "${ME_ENV}" python -m pip install MobileElementFinder

# 2) IntegronFinder env
conda create -y -n "${IF_ENV}" python=3.11 integron_finder

# 3) ISEScan env
conda create -y -n "${IS_ENV}" isescan

# Quick smoke tests
set +e
conda run -n "${ME_ENV}" mefinder find --help >/dev/null 2>&1
ME_OK=$?
conda run -n "${IF_ENV}" integron_finder --version >/dev/null 2>&1
IF_OK=$?
conda run -n "${IS_ENV}" isescan.py --version >/dev/null 2>&1
IS_OK=$?
set -e

echo "Install check:"
echo "  MobileElementFinder: $([[ ${ME_OK} -eq 0 ]] && echo OK || echo FAIL)"
echo "  IntegronFinder:     $([[ ${IF_OK} -eq 0 ]] && echo OK || echo FAIL)"
echo "  ISEScan:            $([[ ${IS_OK} -eq 0 ]] && echo OK || echo FAIL)"

echo
echo "Environments created:"
echo "  ${ME_ENV}"
echo "  ${IF_ENV}"
echo "  ${IS_ENV}"
