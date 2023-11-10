# FLEXR
![logo](img/logo.png)

`FLEXR-GUI` is a Coot1 plugin for [FLEXR](https://github.com/thefischerlab/flexr) - a multiconformer modeling tool.

If you use this software, please cite:

Stachowski, T. R. & Fischer, M.
[FLEXR: automated multi-conformer model building using electron-density map sampling.](https://doi.org/10.1107/S2059798323002498)
2023 Acta Cryst. D79.

## Installation

`FLEXR-GUI` was tested on Intel and M1 Macs running macOS Sonoma.
You will need the following tools installed and accessible in your path:

1. git
2. [Phenix](https://phenix-online.org) (Required to run Ringer)
3. [Homebrew](https://brew.sh)
4. [Coot1](https://github.com/pemsley/coot)

Once these are installed, you need to install these python libraries to the Homebrew python:
1. Install Python libraries (approximate path for M1 Macs):

```
/opt/homebrew/bin/python3.11 -m pip install pandas numpy scipy matplotlib
```

2. Clone the latest release of `FLEXR-GUI`:

```
git clone https://github.com/TheFischerLab/FLEXR-GUI.git
```

3. Move the contents of `FLEXR-GUI` into the Coot1 Python directory:

On M1 Macs:
```
/opt/homebrew/Cellar/coot/<version>/lib/python3.11/site-packages/coot
```

On Intel Macs:
```
/user/local/homebrew/Cellar/coot/<version>/lib/python3.11/site-packages/coot
```

4. Launching Coot1 with the `FLEXR-GUI`:

```
/path/to/bin/coot --script /path/to/flexr_extentions.py
```

5. (optional) A variable for this can also be created in your path with:
```
export coot1='/path/to/bin/coot --script /path/to/flexr_extentions.py'
```

## Usage

1. When you open Coot1 - you should see a button for the `FLEXR-GUI` menu:

![FLEXR-GUI1](img/flexr-gui1.png)


2. FLEXR can be run in two ways - with and without precomputed Ringer measurements.

2a. If you want to run FLEXR, including Ringer, then you need: (1) a PDB, (2) an MTZ, and (3) Phenix installed.
![FLEXR-GUI1](img/flexr-gui2.png)

2b. If you you have precomputed Ringer measurements, then you need (1) a PDB and (2) the path to the Ringer output CSV.
![FLEXR-GUI1](img/flexr-gui3.png)

## Output

1. Running the `FLEXR-GUI` will create a new PDB file with alternative side chain conformations (if any are found).
New alt confs can be quickly assessed with a native Coot1 `altconf GUI` (Draw->Molecule->Residues with Alt Confs...):

![FLEXR-GUI](img/flexr-gui4.png)

2. The typical `FLEXR` output files are saved in the directory from which Coot1 was launched.
More information on that can be found in the original [FLEXR](https://github.com/thefischerlab/flexr) repository.
