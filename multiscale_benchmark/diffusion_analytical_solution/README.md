# Code to compute and plot the analytical solution

This folder contains the implementation of analytical solution and a script to plot it.
Two directories are provided:
- `Decay/` - for case of diffusion with homogeneous uptake
- `CellsAtCenter/` - for case of 12 cells as point sinks at fixed positions 

Each folder provides:
- `AnalyticalSolution.py` - a class representing the analytical solution with helper methods to compute
                            eigenvalues, eigenfunctions, and most importantly amount() and concentration()
                            to compute respectively amount on a given box's volume at a given point in space
                            and concentration at a given point in space.
- `plot_analytical.py` - a script to plot the analytical solution. 
                         This plots average concentration over time, concentration at central point (0,0,0) over time
                         and a 2D plot at stationary state in (x,y) plane crossing at z=0.
- `requirements.txt` - provides the necessary libraries to execute the code

# Usage

Set up your `venv` virtual environment as
```terminal
python3 -m venv path/to/my_env
```

Then activate it by running 
```terminal
source path/to/my_env/bin/activate
```

Once in the virtual environment, install the required librairies as provided by `requirements.txt`
```terminal
python3 -m pip install -r requirements.txt
```

To generate plots, one can simply run with the mentioned options
```terminal
plot_analytical.py [-h] [-p PREFIX] [-D DIFFUSION] [--rate RATE] [-T TEND] [--count-modes COUNT_MODES] [--without-center] [--without-average] [--without-stationary-2D]
```

To display help with detailed options just run
```terminal
plot_analytical.py -h
```
