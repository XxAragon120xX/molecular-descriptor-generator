# 3D Molecular Descriptors Analysis

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![RDKit](https://img.shields.io/badge/RDKit-2022.9+-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Active-brightgreen.svg)

*Comprehensive analysis of molecular descriptors using multiple 3D conformation methods*

---

## Table of Contents

- [Description](#description)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Methodology](#methodology)
- [Project Structure](#project-structure)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Description

This project implements a complete pipeline for 3D molecular descriptors analysis using RDKit and Mordred libraries. The script generates multiple molecular conformations using different geometry algorithms and calculates robust distribution statistics for over 200 molecular descriptors.

### Why is this project useful?

- **Robust Analysis**: Uses 6 different conformation generation methods
- **Robust Statistics**: Performs multiple iterations to obtain statistical distributions
- **Complete Descriptors**: Combines 2D and 3D descriptors for comprehensive analysis
- **Automation**: Fully automated pipeline from SMILES to final descriptors

---

## Features

### 3D Conformation Methods
- **UFF** (Universal Force Field) - Universal force field
- **MMFF** (Merck Molecular Force Field) - Merck molecular force field
- **ETDG** - Distance geometry with experimental torsional preferences
- **KDG** - Knowledge-based distance geometry
- **ETKDGv1** - ETDG + KDG combination version 1
- **ETKDGv2** - ETDG + KDG combination version 2 (improved)

### Types of Calculated Descriptors
- **Surface Descriptors**: PNSA, PPSA, DPSA, FNSA, FPSA, WNSA, WPSA
- **Geometric Descriptors**: Diameter, Radius, Shape index, Petitjean index
- **Gravitational Descriptors**: GRAV, GRAVH, GRAVp, GRAVHp
- **3D-MoRSE Descriptors**: 160 descriptors with different weightings
- **Moments of Inertia**: MOMI-X, MOMI-Y, MOMI-Z, PBF

### Statistical Analysis
- **Multiple Iterations**: 50 repetitions by default for statistical robustness
- **Descriptive Statistics**: Mean, standard deviation, percentiles, min/max
- **Reproducibility**: Controllable random seed

---

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)

### Dependencies Installation

```bash
# Clone the repository
git clone https://github.com/your-user/molecular-descriptors-analysis.git
cd molecular-descriptors-analysis

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Main Dependencies

```bash
pip install rdkit-pypi>=2022.9.1
pip install pandas>=1.5.0
pip install numpy>=1.21.0
pip install mordred>=1.2.0
pip install matplotlib>=3.5.0
pip install openpyxl>=3.0.9
pip install py3Dmol>=1.8.0
```

---

## Usage

### Basic Usage

```python
python molecular_analysis.py
```

### Input File Structure

Your Excel file should have the following structure:

| SMILES | Name | Property1 | ... |
|--------|------|-----------|-----|
| CCO | Ethanol | 0.789 | ... |
| c1ccccc1 | Benzene | 0.876 | ... |
| CC(=O)O | Acetic acid | 1.049 | ... |

**File requirements:**
- Format: `.xlsx` (Excel)
- File name: `DATA.xlsx` (or modify in code)
- Required column: `SMILES` (with valid SMILES notation)

### Customization

```python
# Modify number of iterations
number_iterations = 100  # Default: 50

# Change input file
df_original = pd.DataFrame(pd.read_excel('my_dataset.xlsx'))

# Select specific descriptors
descriptor_names_3d = ["PNSA1", "PPSA1", "GeomDiameter"]
```

---

## Methodology

### Processing Pipeline

```
SMILES File → RDKit Molecules → 3D Conformations → Descriptors Calculation → Statistical Analysis → Data Cleaning → Final Export
```

### Detailed Flow

1. **Data Loading**: Reads SMILES structures from Excel
2. **Conformation Generation**: Creates 3D geometries with 6 different methods
3. **Iterative Calculation**: Repeats the process N times for statistical analysis
4. **Statistical Analysis**: Calculates descriptive statistics per descriptor
5. **Data Combination**: Integrates 2D descriptors and 3D statistics
6. **Cleaning**: Removes non-numeric values and handles missing data
7. **Export**: Saves results in Excel format

---

## Project Structure

```
molecular-descriptors-analysis/
├── README.md                        # This file
├── molecular_analysis.py            # Main script
├── requirements.txt                 # Project dependencies
├── data/                           # Input data
│   ├── DATA.xlsx                   # Example file
│   └── example_dataset.xlsx        # Test dataset
├── results/                        # Generated results
│   └── complete_molecular_descriptors_analysis.xlsx
├── docs/                           # Additional documentation
│   ├── methodology.md              # Detailed methodology
│   └── descriptors_reference.md    # Descriptors reference
└── tests/                          # Unit tests
    └── test_analysis.py            # Main script tests
```

---

## Results

### Output File

The script generates an Excel file with the following characteristics:

- **Original Columns**: Maintains all columns from original dataset
- **2D Descriptors**: ~1800 two-dimensional descriptors
- **3D Statistics**: Distribution statistics for each 3D descriptor:
  - Mean (`mean`)
  - Standard deviation (`std`)
  - 25th percentile (`25%`)
  - Median (`50%`)
  - 75th percentile (`75%`)
  - Minimum (`min`)
  - Maximum (`max`)

### Example of Generated Columns

```
Original: SMILES, Name, Activity, ...
2D: nAcid, nBase, MolWt, TPSA, ...
3D: PNSA1_mean_ETKDGv2, PNSA1_std_ETKDGv2, GeomDiameter_mean_UFF, ...
```

### Quality Metrics

- **Completeness**: >95% of molecules processed successfully
- **Speed**: ~1-2 molecules/second (hardware dependent)
- **Precision**: >99% reproducibility with same random seed

---

## Contributing

Contributions are welcome! Here's how you can help:

### Ways to Contribute

1. **Report Bugs**
   - Use the issue tracker
   - Include detailed error information
   - Provide reproducible example

2. **Suggest Improvements**
   - New conformation methods
   - Additional descriptors
   - Performance optimizations

3. **Submit Pull Requests**
   - Fork the repository
   - Create a branch for your feature
   - Add tests if necessary
   - Submit the pull request

### Development Guide

```bash
# Fork and clone
git clone https://github.com/your-user/molecular-descriptors-analysis.git

# Create branch for new feature
git checkout -b feature/new-feature

# Make changes and commit
git add .
git commit -m "Add new feature"

# Push and create pull request
git push origin feature/new-feature
```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

```
MIT License

Copyright (c) 2024 Your Name

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software...
```

---

## Contact

### Author
**Your Name**
Daniel Aragon Giraldo
- Email: danielaragon120@gmail.com
- GitHub: [@XxAragon120xX](https://github.com/XxAragon120xX)
- LinkedIn: [Daniel Aragon Giraldo](https://www.linkedin.com/in/daniel-aragon-giraldo-2bb415344/)

### Support
- **Report Bug**: Create Issue
- **Suggest Feature**: Create Issue
- **Questions**: Discussions

---

## Acknowledgments

This project uses the following open-source libraries:

- **[RDKit](https://www.rdkit.org/)** - Cheminformatics toolkit
- **[Mordred](https://github.com/mordred-descriptor/mordred)** - Molecular descriptors library
- **[Pandas](https://pandas.pydata.org/)** - Data analysis in Python
- **[NumPy](https://numpy.org/)** - Scientific computing

---

## Project Statistics

![GitHub repo size](https://img.shields.io/github/repo-size/your-user/molecular-descriptors-analysis)
![GitHub last commit](https://img.shields.io/github/last-commit/your-user/molecular-descriptors-analysis)
![GitHub issues](https://img.shields.io/github/issues/your-user/molecular-descriptors-analysis)
![GitHub pull requests](https://img.shields.io/github/issues-pr/your-user/molecular-descriptors-analysis)

---

**If this project has been useful to you, don't forget to give it a star**

[Back to top](#3d-molecular-descriptors-analysis)
