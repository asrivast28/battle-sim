Battle Simulator using Floor Field Models
=========================================
(CSE6730 Modeling and Simulation Project 2)


Building
--------

Our project consists of two parts: a python frontend and a C++ backend. The
python frontend enables fancy, interactive visualizations of the simulated
battles. The backend implements the actual floor field model. The parallelized
floor field model is implemented only in C++, not in the python frontend.

### Parallelized Floor Field

Build requirements:
- MPI (MPI-2 capable MPI implementation, such as OpenMPI version 1.6.5)
- cmake

Building:
```sh
mkdir build
cd build
cmake ../
make
```

### Python frontend (visualizations)

Build requirements:
- python 2.7
- pygames
- moviepy
- swig
- gizeh

Building:
```sh
cd driver
make clean && make
```

Running the code
----------------

### Visualizations

We provide two different visualization engines. One interactive one, and one
which can render movie and `gif` files. The following commands
have to be executed inside the `driver` folder:

Rendering videos/gifs with `./simulator`:
```
usage: simulator [-h] [--seed SEED] [--iter COUNT] [--fps FPS]
                 [--render CHOICE] [--ext EXT]

optional arguments:
  -h, --help       show this help message and exit
  --seed SEED      Random seed value for running the simulation.
  --iter COUNT     Number of iterations. One iteration includes move and kill.
  --fps FPS        Frames per second.
  --render CHOICE  Choice of rendering used by the simulator {pixel, icon}.
  --ext EXT        File type to be generated {mp4, gif}.
```

To run the interactive visualization (game), run:
```
./game
```

then you can draw the two different armies and set the targets. Drawing is
only possible in *edit* mode, whereas setting the target is always possible.
To draw or set the target, simply press the left mouse button into the window.
The key bindings are as follows:
- **b**: Draw blue army
- **g**: Draw green army
- **h**: Set green armies target
- **n**: Set blue armies target
- **s**: Draw static field (walls and such)
- **SPACE**: Toggle between *edit* and *play* mode.

