#!/usr/bin/env python3
import sys
from types import SimpleNamespace
import argparse
import numpy as np


def parse_arguments(radii):
    p = argparse.ArgumentParser()
    p.add_argument('--num', action='store_true',  help='add atom numbers')
    p.add_argument('--canvas-size', default=80.0, type=float,  help='basic canvas size')
    p.add_argument('-wa', '--atom-border', default=5.0, type=float,  help='atom border width (default 5.0)')
    p.add_argument('-wb', '--bond-width', default=5.0, type=float,  help='bond line width (default 5.0)')
    p.add_argument('-db', '--bond-distance', default=0.05, type=float,  help='line distance in a multiple bond (default 0.05)')
    p.add_argument('-rs', '--atom_size', default=1.0, type=float,  help='scaling factor for atom radii')
    for i in range(len(radii)):
        p.add_argument(f'-r{i}', default=None, type=float, help=('-r{q} sets radius for element {q} in Å' if i==6 else argparse.SUPPRESS))
    p.add_argument('-g',  '--gradient', action='store_true', help='fill atoms with radial gradients')
    p.add_argument('-fs', '--font-size', default=16, type=int, help='font size (default 16)')
    p.add_argument('-fn', '--font-name', default='monospace', type=str, help='font name (default monospace)')
    args = p.parse_args()

    for i in range(len(radii)):
        r = eval(f'args.r{i}')
        if r is not None:
            radii[i] = r
    radii *= args.atom_size * 0.2

    bond_style = SimpleNamespace(
            color = 'black',
            width = args.bond_width,
            distance = args.bond_distance,
            )
    atom_style = SimpleNamespace(
            stroke_color = 'black',
            stroke = args.atom_border,
            )
    text_style = SimpleNamespace(
            font = args.font_name,
            size = args.font_size,
            stroke = 1.,
            fill_color = '#FFFFFF',
            stroke_color = '#000000',
            stroke_opacity = 1.0,
            fill_opacity = 1.0,
            weight = 'bold',
            style = 'normal',
            )
    par = SimpleNamespace(text=text_style,
                          atom=atom_style,
                          bond=bond_style,
                          num=args.num,
                          canvas_size=args.canvas_size,
                          grad=args.gradient)
    return par


def print_svg(Q, R, bonds, labels, radii, colors, colors_ini, colors_fin, par):

    radmax = max(radii[Q])
    rmin = np.array((min(R[:,0]), min(R[:,1])))
    rmax = np.array((max(R[:,0]), max(R[:,1])))
    center = np.array((rmin[0]-radmax, rmax[1]+radmax, 0.0))
    rsize = rmax - rmin + 2*radmax
    z_order = np.argsort(R[:,2])

    a = par.canvas_size
    afactor = np.array((a, -a, 1.0))
    rcanv = (a*rsize).astype(int)

    print("<svg xmlns=\"http://www.w3.org/2000/svg\" "
          "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
          f"width=\"{rcanv[0]}\" height=\"{rcanv[1]}\">\n")

    print("  <defs>")
    if par.grad:
        for q in sorted(set(Q)):
                print(f'    <g id="atom{q}"> '
                      f'<radialGradient id="radGrad{q}" cx="0.5" cy="0.5" fx="0.33" fy="0.33" r="0.66"> '
                      f'<stop offset="0%" stop-color="#{colors_ini[q]:06x}" /> '
                      f'<stop offset="100%" stop-color="#{colors_fin[q]:06x} "/> '
                      f'</radialGradient> '
                      f'<circle cx="0" cy="0" r="{abs(radii[q])*a}" '
                      f'fill="url(#radGrad{q})"/>'
                      '</g>')
    else:
        for q in sorted(set(Q)):
                print(f'    <g id="atom{q}"> '
                      f'<circle cx="0" cy="0" r="{abs(radii[q])*a}" '
                      f'fill="#{colors[q]:06x}" stroke="{par.atom.stroke_color}" stroke-width="{par.atom.stroke*a/132.317536}"/> '
                      '</g>')
    print("  </defs>\n")

    if par.num or labels:
        print("  <style>\n"
            "  .atnum {\n"
            f"   font-size:{par.text.size*a/80.0}px;font-style:{par.text.style}; font-weight:{par.text.weight};\n"
            f"   font-family:'{par.text.font}';\n"
            f"   fill:{par.text.fill_color}; fill-opacity:{par.text.fill_opacity}; \n"
            f"   stroke:{par.text.stroke_color};stroke-width:{par.text.stroke*a/80.0}; stroke-opacity:{par.text.stroke_opacity};\n"
            "    text-align:start; writing-mode:lr-tb; text-anchor:start;\n"
            "  }\n"
            "  </style>\n")


    for i, I in enumerate(z_order):

        ri = afactor * (R[I]-center)
        print(f'  <use x="{ri[0]}" y="{ri[1]}" xlink:href="#atom{Q[I]}"/>')

        for J in z_order[i+1:]:

            nb = bonds[I,J]
            if nb==0:
                continue
            dash = '' if nb > 0 else ' stroke-dasharray="10,5"'

            rij = R[J]-R[I]
            rij *= 0.666 / np.linalg.norm(rij)

            ri = afactor * (R[I]-center + rij*abs(radii[Q[I]]))
            rj = afactor * (R[J]-center - rij*abs(radii[Q[J]]))

            for ib in range(-abs(nb)+1, abs(nb), 2):
                x1 = ri[0] + rij[1]*ib*a*par.bond.distance
                y1 = ri[1] + rij[0]*ib*a*par.bond.distance
                x2 = rj[0] + rij[1]*ib*a*par.bond.distance
                y2 = rj[1] + rij[0]*ib*a*par.bond.distance
                print(f'    <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{par.bond.color}" stroke-width="{par.bond.width*a/132.317536}"{dash}/>')

    if par.num:
        print()
        for i, ri in enumerate(R):
            ri = afactor * (ri-center)
            print(f'  <text x="{ri[0]}" y="{ri[1]}" class="atnum">{i+1}</text>')
    if labels:
        print()
        for i, label in labels:
            ri = R[i]
            ri = afactor * (ri-center)
            print(f'  <text x="{ri[0]}" y="{ri[1]}" class="atnum">{label}</text>')

    print('</svg>')

def mol_input():
    data = [*map(lambda x: x.strip().split(), sys.stdin.readlines())]

    Q = []
    R = []
    for atom in filter(lambda x: x[0]=='atom', data):
        Q.append(int(atom[1]))
        R.append([*map(float, atom[2:5])])

    bonds = np.zeros((len(Q), len(Q)), dtype=int)
    for bond in filter(lambda x: x[0]=='bond', data):
        i,j = map(lambda x: int(x)-1, bond[1:3])
        order = 1 if len(bond)<=3 else int(bond[3])
        bonds[i,j] = bonds[j,i] = order

    labels = []
    for label in filter(lambda x: x[0]=='label', data):
        i,text = int(label[1])-1, label[2]
        labels.append((i, text))
    if len(labels)==0:
        labels = None

    return np.array(Q), np.array(R), bonds, labels


def atom_parameters():
    radii = [ 0.0,
      0.35, 0.20,
      1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.40,
      1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 0.90,
      2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35,
      1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.00,
      2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40,
      1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.30,
      2.35, 1.98, 1.69, 1.65, 1.65, 1.64, 1.65, 1.66, 1.85,
                  1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.70,
                  1.56, 1.44, 1.34, 1.30, 1.28, 1.26, 1.26, 1.29,
      1.34, 1.44, 1.55, 1.54, 1.52, 1.53, 1.53, 1.50,
      2.70, 2.23, 1.87, 1.78, 1.61, 1.40, 1.40, 1.40, 1.40,
                  1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40,
                  1.40, 1.40
           ] + [1.40]*23
    colors = [ 0x000000,
      0x808080, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xC0C0C0, 0x0000FF, 0xFF0000, 0xFF00FF, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0x090909, 0x0F0B0B, 0xF0F000, 0x00F000, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
            ] + [0xA0A0A0]*23
    colors_ini = [ 0x000000,
      0xd0d0d0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xCC22CC, 0x808080, 0x2222FF, 0xFF2222, 0xFF22FF, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0x090909, 0x0F0B0B, 0xFFFF22, 0x22FF22, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0,
                          0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0xA0A0A0,
            ] + [0xA0A0A0]*23
    colors_fin = [ 0x000000,
      0x606060, 0x303030,
      0x303030, 0x303030, 0x552255, 0x202020, 0x000088, 0x880000, 0x880088, 0x303030,
      0x303030, 0x303030, 0x303030, 0x090909, 0x0F0B0B, 0x888800, 0x008800, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
                          0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
                          0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
      0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
                          0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
                          0x303030, 0x303030,
                          0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030, 0x303030,
            ] + [0x303030]*23
    return np.array(radii), np.array(colors), np.array(colors_ini), np.array(colors_fin)


if __name__=='__main__':
    radii, colors, colors_ini, colors_fin = atom_parameters()
    parameters = parse_arguments(radii)
    Q, R, bonds, labels = mol_input()
    print_svg(Q, R, bonds, labels, radii, colors, colors_ini, colors_fin, parameters)

