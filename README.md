# Quantum Assembly Theory (QAT) Framework ⚛️

[![ChemRxiv](https://img.shields.io/badge/ChemRxiv-Preprint-29b6f6?style=for-the-badge&logo=arxiv&logoColor=white)](#) 
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=for-the-badge&logo=python&logoColor=white)](#)

> **Official Supplementary Information Repository** > **Author:** Akshat Sharma (IISER Kolkata)  
> **Preprint:** *Quantum assembly theory defines molecular complexity from electron density* (ChemRxiv 2026) 

## 📋 Overview
This repository contains the full computational pipeline used to calculate the Density-Based Assembly Index (Aρ). The workflow integrates semi-empirical geometry optimization, DFT single-point calculations, QTAIM topological analysis, and custom Python strain calculators to evaluate molecular complexity from fundamental electron density descriptors.

---

## 🚀 The Aρ Computational Pipeline

The framework is executed in a strict five-step protocol:

### Step 1: Optimization & Wavefunction Generation
Initial structural coordinates (`.sdf` or `.xyz`) are refined before wavefunction generation.
1. **Geometry Optimization:** Performed using **xTB**.
2. **Single Point Calculation:** DFT evaluation using the **r2SCAN-3c** functional.
3. **Conversion:** The final electronic structure output is converted into a standard `molden.input` file for topological analysis.

### Step 2: Topological Analysis (Multiwfn)
Critical points are located and evaluated using **Multiwfn**. 
Load the `molden.input` file into Multiwfn and execute the following exact command sequence for automated batch processing:
```text
2, 2, 3, 4, 5, 8, 7, 0, 0
