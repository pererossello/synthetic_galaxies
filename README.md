# Galactic Particle Orbit Visualizer

Generates videos of particle orbits in various galactic potential models. 

<!-- ![Video GIF Placeholder](figures/NFW_k/video_20.gif)
![Alt text](figures/NFW_k/video_20.gif "NFW") -->

![One million particles under a Navarro-Frenk-White potential](NFW_example.gif)


## Structure
- `code/gutils.py`: 
  - Functions to generate density distributions.
  - Sample particle positions.
  - Compute initial velocities.

- `code/plot_utils.py`: 
  - Utility functions for plotting the orbits.

- `code/compute_and_plot.ipynb`: 
  - Main Jupyter notebook for generating frames at each time step of the integrated orbits.

## Dependencies
This project requires the following Python packages:
- [Gala](https://github.com/adrn/gala): A Python library for galactic dynamics.
- [Astropy](https://www.astropy.org/): A community-developed core Python package for Astronomy.
- [SymPy](https://www.sympy.org/): A Python library for symbolic mathematics.

## License
This project is licensed under the [MIT License](LICENSE.md).
