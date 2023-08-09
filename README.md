# xyz2svg

Lightweight script to make vector images of molecules

## Requirements
- `python>=3.6`
- [briling/v](https://github.com/briling/v) or something capable of printing atom positions and connectivity


## Basic usage

```
v mol/caffeine.xyz | ./mol2svg.py > caffeine.svg
```
Open an xyz file with `v`, rotate the molecule as you like and press `p`, then `x` to exit

![vector balls-and-sticks image of caffeine](figures/caffeine.svg)

## Input file modification
For a water dimer
```
v mol/kj-H2O_1--H3O+.xyz > input/kj-H2O_1--H3O+.inp
```
prints
```
atom   8   -1.1697501   -0.1527130    0.0202511
atom   1   -1.8639331    0.4812310   -0.1588469
atom   1   -1.6335311   -0.9640680    0.2263131
atom   8    1.1960949    0.3707430   -0.0773079
atom   1    0.2351029    0.2036530   -0.0496909
atom   1    1.6046709    0.2360230    0.7985351
atom   1    1.6313459   -0.1748690   -0.7592539
bond   1   2
bond   1   3
bond   4   5
bond   4   6
bond   4   7
```
This file contains atomic numbers and coordinates in Å and list of bonds (base 1).
We can add a hydrogen bond between atoms 1–5:
```
bond 1 5 -1
```
(see [input/kj-H2O_1--H3O+.inp](input/kj-H2O_1--H3O+.inp)). Then 
```
./mol2svg.py < input/kj-H2O_1--H3O+.inp > figures/kj-H2O_1--H3O+.svg
```
gives<br>
![vector balls-and-sticks image of a water dimer](figures/kj-H2O_1--H3O+.svg)

One can also add (positive or negative) bond orders ([input/caffeine_v2.inp](input/caffeine_v2.inp)):
![vector balls-and-sticks image of caffeine with double bonds](figures/caffeine_v2.svg)

## Drawing options

Run 
```
./mol2svg.py --help
```
to get help:
```
  --num                                                 add atom numbers
  --canvas-size CANVAS_SIZE                             basic canvas size (default 80)
  -wa ATOM_BORDER, --atom-border ATOM_BORDER            atom border width (default 5.0)
  -wb BOND_WIDTH, --bond-width BOND_WIDTH               bond width (default 5.0)
  -db BOND_DISTANCE, --bond-distance BOND_DISTANCE      line distance in multiple bonds (default 0.05)
  -rs ATOM_SIZE, --atom_size ATOM_SIZE                  scaling factor for atom radii
  -r%d ATOM_%d_RADIUS                                   sets basic radius for a specific element (in Å)
  -g, --gradient                                        fill atoms with radial gradients (pseudo-3D mode)
```
The same caffeine with options:
```
./mol2svg.py -rs 3.5 -wa 10 -wb 10 -db 0.1 < input/caffeine_v2.inp > figures/caffeine_v3.svg
```
![vector balls-and-sticks image of caffeine with double bonds and looks nice](figures/caffeine_v3.svg)


@iribirii's pseudo-3D mode:
```
./mol2svg.py -rs 3.5 -wb 10 -g < input/caffeine.inp > figures/caffeine_v4.svg
```
![vector balls-and-sticks image of caffeine with a pseudo 3D effect](figures/caffeine_v4.svg)
