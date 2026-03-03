# Project_2026_I

List of participants:

- JONATHAN BARRIENTOS ANDRADE (J-dot-Barrientos)
- ADRIÀ BRÚ I CORTÉS (cooligula) Main Leader
- MARC FREIXER REALP (Fresco4)
- ADRIAN LLAMAS JARAMILLO (AdrianLLJ)

List of tasks: Responsible of task
- Initialization: JONATHAN BARRIENTOS ANDRADE
- Energy: ADRIAN LLAMAS JARAMILLO
- Monte Carlo Update: MARC FREIXER REALP
- Post processing and statistics: ADRIÀ BRÚ I CORTÉS

## Project structure
```
mc-polymer-simulation/
├── Makefile             # Main Makefile
├── README.md            
├── .gitignore           
├── src/                 # Directory containing all source code.
│   ├── sequential/      # Initial sequential code for the March 12th deliverable.
│   └── parallel/        # The MPI version supporting replica, parallel energy, and mixed parallelism.
├── scripts/             # Python scripts for post-processing and plotting.
├── results/             
│   ├── data/            # Energy, torsion angle, etc.
│   ├── plots/           # Plots of obtained data.
│   └── trajectory/      # Monte Carlo trajectory.
└── tests/               # Differents tests.
```