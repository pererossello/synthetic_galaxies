# Galactic Particle Orbit Visualizer

This project generates videos of particle orbits in various galactic potential models. It leverages the `gala` Python package to integrate the orbits and create dynamic visualizations.

## Demo Video
![Video GIF Placeholder](path/to/video.gif)

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

## Citation
If you use this software in your research, please cite it as follows:
