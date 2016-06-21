Qahwa
====

command line tools for doing Umbrella Sampling & WHAM.

The main purpose of Qahwa is calculating Potential of Mean Force(PMF)
based on dcd files including umbrella sampling trajectory.

The sample result is shown below.
![umbrella](http://10.1.1.222/gitlab/niina/Qahwa/raw/master/sample/umbrella_biased.png)

Fig.1 : The umbrella sampling trajectory

![wham](http://10.1.1.222/gitlab/niina/Qahwa/raw/master/sample/wham_result.png)

Fig.2 : The result of wham

## Build

clone this repository on your machine.

    $ git clone git@moca:niina/Qahwa.git

At first, you should initialize submodules.

    $ git submodule update --init --recursive

then use CMake.

    $ cd build
    $ cmake ..
    $ make

## Usage

### Reaction Coordinate and PerturbingPotential

To use Qahwa, you must write very short c++ code to calculate your
Reaction Coordinate and Perturbing Potential from SnapShot.

the class 'SnapShot' is an alias of std::vector<Vector3d>, and Vector3d is
normal 3-dimentional vector defined in AX library.
std::vector is a kind of dynamic array.

You can access the position of i-th particle like this way.

    // qahwa code
    SnapShot snapshot;
    Vector3d position = snapshot.at(i);

__NOTE__: Because the c++ array is 0-based index,
          if you want to access the position of particle that the imp is n,
          you should access at(n-1) element of snapshot.

And the AX library provides you some useful functions for vectors. like

    // qahwa c++ code
    Vector3d position(0.0, 1e0, 1.5e2);
    Vector3d another_position = snapshot.at(i);
    double   l  = length(position);
    double   l2 = len_square(position);
    double   dot = dot_prod(position, another_position);
    Vector3d cross = cross_prod(position, another_position);

As default, the harmonic anchor is implemented as PerturbingPotential and
distance between first 2 particle is implemented as ReactionCoordinate.

See __./src/UserDefinedFunction.hpp__ and __./src/UserDefinedFunction.cpp__

### Command line interface

Qahwa has several modes.

| mode             | input file | description                                                         |
|:-----------------|:-----------|:--------------------------------------------------------------------|
| --make-histogram | dcd file   | make histogram from one dcd file(all the values are set as default) |
| --make-histogram | toml file  | make histogram from several dcd file                                |
| --make-unbiased  | toml file  | make unbiased histogram of each trajectories (not connected)        |
| --make-pmf       | toml file  | make potential of mean force from probability density function      |
| --wham           | toml file  | run wham and output weighted unbiased histogram of the trajs        |

To use qahwa, run this way.

    $ qahwa --wham input.toml

### input file

Input file requires toml format.

The typical example is described below.

    [histogram]
    dcdfiles   = ["traj1.dcd", "traj2.dcd", "traj3.dcd"]
    output     = "histogram"

    [unbiased]
    dcdfiles   = ["traj1.dcd", "traj2.dcd", "traj3.dcd"]
    parameters = [
        [1, 1.0, 5.0, 0.0, 0.0],
        [1, 1.0, 5.5, 0.0, 0.0],
        [1, 1.0, 6.0, 0.0, 0.0]
    ]
    output     = "unbiased"

    [wham]
    dcdfiles   = ["traj1.dcd", "traj2.dcd", "traj3.dcd"] # required
    parameters = [
        [1, 1.0, 5.0, 0.0, 0.0],
        [1, 1.0, 5.5, 0.0, 0.0],
        [1, 1.0, 6.0, 0.0, 0.0]
    ] # required
    bins        = 200    # default is 100.
    temperature = 300.0  # default is 300.0
    output      = "wham" # required

"dcdfiles" is the list of "./path/to/filename"s.

"parameters" is the list of parameters, and
a parameter is described as a list of floating-point number.
As default, the first element is index of anchored particle
(though it is float...it is a problem)
and the second element is the coefficient of anchor.
And the last 3 element is anchor position(3d-vector).
"output" is a name of the file to output.

Qahwa reads only the table that is related to the mode. For example, if you run
Qahwa with option --make-unbiased, [histogram] and [wham] tables are ignored.

Because of my not implementing toml array\_of\_table feature,
the input format is a bit ugly and uneasy to read.
After I implementing more powerful toml parser,
the format will be changed as more elegant.
for example, like this.

    [wham]
    dcdfiles = ["traj1.dcd", "traj2.dcd", ...]

    [[wham.parameters]]
    index = 1
    coefficient = 1.0
    position = [5.0, 0.0, 0.0]

    [[wham.parameters]]
    index = 1
    coefficient = 1.0
    position = [6.0, 0.0, 0.0]

    ...

## Licensing terms

This project is licensed under the terms of the MIT License.
See LICENSE for the project license.

- Copyright (c) 2016- Toru Niina

All rights reserved.
