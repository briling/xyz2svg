#!/usr/bin/env python3
import sys
from types import SimpleNamespace
import argparse
import numpy as np


def parse_arguments(radii):
    p = argparse.ArgumentParser()
    p.add_argument('--num', action='store_true',  help='add atom numbers')
    p.add_argument('--elements', action='store_true',  help='add element symbols')
    p.add_argument('--canvas-size', default=80.0, type=float,  help='basic canvas size')
    p.add_argument('-wa', '--atom-border', default=5.0, type=float,  help='atom border width (default 5.0)')
    p.add_argument('-wb', '--bond-width', default=5.0, type=float,  help='bond line width (default 5.0)')
    p.add_argument('-db', '--bond-distance', default=0.05, type=float,  help='line distance in a multiple bond (default 0.05)')
    p.add_argument('-rs', '--atom_size', default=1.0, type=float,  help='scaling factor for atom radii')
    for i in range(len(radii)):
        p.add_argument(f'-r{i}', default=None, type=float, help=('-r{q} sets radius for element {q} in Å' if i==6 else argparse.SUPPRESS))
    p.add_argument('-g',  '--gradient', action='store_true', help='fill atoms with radial gradients')
    p.add_argument('--fog', action='store_true', help='enable fog for depth perspective')
    p.add_argument('--fog-strength', default='0.8', type=float, help='fog strength (default 0.8, between 0.0 and 1.0)')
    p.add_argument('--light-hydrogen', action='store_true', help='use a lighter color for H')
    p.add_argument('--no-H', action='store_true', help='hide hydrogens bonded to carbon')
    p.add_argument('-fs', '--font-size', default=24, type=int, help='font size (default 24)')
    p.add_argument('-fn', '--font-name', default='monospace', type=str, help='font name (default monospace)')
    p.add_argument('--bond-color', default='#000000', type=str,  help='bond line color (default black - hex)')
    p.add_argument('--atom-stroke-color', default='#000000', type=str,  help='atom stroke color (default black - hex)')
    p.add_argument('--text-stroke-color', default='#FFFFFF', type=str,  help='text stroke color (default white - hex)')
    p.add_argument('--text-color', default='#000000', type=str,  help='text fill color (default black - hex)')
    p.add_argument('--text-weight', default='bold', type=str,  help='text weight (default bold)')
    p.add_argument('--text-style', default='normal', type=str,  help='text style (default normal)')
    p.add_argument('--text-stroke-width', default=8, type=int,  help='text stroke width (default 8)')
    p.add_argument('--value-gradient', nargs=2, default=['#000000', '#FF0000'], type=str, help='starting and finishing colors for value gradient (default ["#000000", "#FF0000"]')
    p.add_argument('--value-radius', default=0.2, type=float, help='radius of value gradient circles (default 0.2 Å)')
    args = p.parse_args()

    if args.num and args.elements:
        raise RuntimeError

    for i in range(len(radii)):
        r = eval(f'args.r{i}')
        if r is not None:
            radii[i] = r
    radii *= args.atom_size * 0.2

    bond_style = SimpleNamespace(
            color = args.bond_color,
            width = args.bond_width,
            distance = args.bond_distance,
            )
    atom_style = SimpleNamespace(
            stroke_color = args.atom_stroke_color,
            stroke = args.atom_border,
            )
    text_style = SimpleNamespace(
            font = args.font_name,
            size = args.font_size,
            stroke = args.text_stroke_width,
            fill_color = args.text_color,
            stroke_color = args.text_stroke_color,
            stroke_opacity = 1.0,
            fill_opacity = 1.0,
            weight = args.text_weight,
            style = args.text_style,
            paint_order = 'stroke fill markers',
            )
    val_grad = SimpleNamespace(
            gcolors = args.value_gradient,
            rvalue = args.value_radius * args.atom_size * 0.2,
            )
    fog_style = SimpleNamespace(
            fog = args.fog,
            strength = args.fog_strength,
            )

    if not (0<fog_style.strength<1):
        print('warning: using default fog strength value')
        fog_style.strength = 0.8

    par = SimpleNamespace(text=text_style,
                          atom=atom_style,
                          bond=bond_style,
                          num=args.num,
                          elements=args.elements,
                          canvas_size=args.canvas_size,
                          grad=args.gradient,
                          val_grad=val_grad,
                          fog=fog_style,
                          light_hydrogen=args.light_hydrogen,
                          no_H=args.no_H
                          )


    return par

def blend_fog(hex_colour, fog_rgb, strength):
    """
    Blend the given colour with the fog colour based on the strength.
    """
    strength = strength**2
    rgb = np.array([int(hex_colour[i:i+2], 16) for i in (1,3,5)])
    blended = ( 1 - strength ) * rgb + strength * fog_rgb
    return "#{:02x}{:02x}{:02x}".format(*blended.astype(int))


def linear_grad(gcolors, value):
    gcolors = np.array([[int(h[i:i+2], 16) for i in (1, 3, 5)] for h in gcolors])
    rgb = gcolors[1]*value + gcolors[0]*(1.0-value)
    rgb = np.round(rgb).astype(int)
    return f'{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}'

def should_skip_atom(atom_idx, Q, bonds, no_H):
    """Check if atom should be hidden (H bonded to C when no_H flag is set)."""
    if not no_H:
        return False
    if abs(Q[atom_idx]) != 1:  # Not hydrogen
        return False
    # Count bonds and check if all are to carbon
    bonded_atoms = []
    for j, bond_order in enumerate(bonds[atom_idx]):
        if bond_order != 0:
            bonded_atoms.append(j)
    
    # Only hide if H has exactly one bond and it's to carbon
    if len(bonded_atoms) == 1 and abs(Q[bonded_atoms[0]]) == 6:
        return True
    return False

def print_svg(Q, R, bonds, labels, values, radii, colors, colors_ini, colors_fin, elements, par):

    radmax = max(radii[Q])
    rmin = np.array((min(R[:,0]), min(R[:,1])))
    rmax = np.array((max(R[:,0]), max(R[:,1])))
    center = np.array((rmin[0]-radmax, rmax[1]+radmax, 0.0))
    rsize = rmax - rmin + 2*radmax
    z_order = np.argsort(R[:,2])

    if par.fog.fog:
        zrange = max(R[:,2].max() - R[:,2].min(), 1e-6)
        fog_factors = par.fog.strength * (R[:,2] - R[:,2].max()) / zrange
        fog_rgb = np.array([255, 255, 255])
    else:
        fog_factors = np.zeros(len(R))

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


    if len(values)>0:
        print(f'    <g id="values"> '
              f'<circle cx="0" cy="0" r="{par.val_grad.rvalue*a}" '
              f'stroke="{par.atom.stroke_color}" stroke-width="{par.atom.stroke*a/132.317536}"/> '
              '</g>')
        pass

    print("  </defs>\n")

    if par.num or par.elements or labels:
        print("  <style>\n"
            "  .atnum {\n"
            f"   font-size:{par.text.size*a/80.0}px;font-style:{par.text.style}; font-weight:{par.text.weight};\n"
            f"   font-family:'{par.text.font}'; paint-order:{par.text.paint_order};\n"
            f"   fill:{par.text.fill_color}; fill-opacity:{par.text.fill_opacity}; \n"
            f"   stroke:{par.text.stroke_color};stroke-width:{par.text.stroke*a/80.0}; stroke-opacity:{par.text.stroke_opacity};\n"
            "    text-align:start; writing-mode:lr-tb; text-anchor:start;\n"
            "  }\n"
            "  </style>\n")


    for i, I in enumerate(z_order):

        if should_skip_atom(I, Q, bonds, par.no_H):
            continue

        ri = afactor * (R[I]-center)
        if par.fog.fog and not par.grad:
            fog_color = blend_fog(f'#{colors[Q[I]]:06x}', fog_rgb, fog_factors[I])
            stroke_color = blend_fog(par.atom.stroke_color, fog_rgb, fog_factors[I])
            print(f'  <circle cx="{ri[0]}" cy="{ri[1]}" '
              f'r="{abs(radii[Q[I]])*a}" '
              f'fill="{fog_color}" '
              f'stroke="{stroke_color}" '
              f'stroke-width="{par.atom.stroke*a/132.317536}"/>')
        else:
            print(f'  <use x="{ri[0]}" y="{ri[1]}" xlink:href="#atom{Q[I]}"/>')

        for J in z_order[i+1:]:
            
            if should_skip_atom(I, Q, bonds, par.no_H) or should_skip_atom(J, Q, bonds, par.no_H):
                continue

            nb = bonds[I,J]
            if nb==0:
                continue
            dash = '' if nb > 0 else ' stroke-dasharray="10,5"'

            rij = R[J]-R[I]
            rij *= 0.666 / np.linalg.norm(rij)

            ri = afactor * (R[I]-center + rij*abs(radii[Q[I]]))
            rj = afactor * (R[J]-center - rij*abs(radii[Q[J]]))
            if par.fog.fog and not par.grad:
                stroke_color = blend_fog(par.bond.color, fog_rgb, (fog_factors[I] + fog_factors[J]) / 2)
            else:
                stroke_color = par.bond.color
            #print('xxx', *R[I], *R[J])
            for ib in range(-abs(nb)+1, abs(nb), 2):
                x1 = ri[0] + rij[1]*ib*a*par.bond.distance
                y1 = ri[1] + rij[0]*ib*a*par.bond.distance
                x2 = rj[0] + rij[1]*ib*a*par.bond.distance
                y2 = rj[1] + rij[0]*ib*a*par.bond.distance
                print(f'    <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke_color}" stroke-width="{par.bond.width*a/132.317536}"{dash}/>')

    if par.num:
        print()
        for i, ri in enumerate(R):
            ri = afactor * (ri-center)
            print(f'  <text x="{ri[0]}" y="{ri[1]}" class="atnum">{i+1}</text>')
    if par.elements:
        print()
        for qi, ri in zip(Q, R):
            ri = afactor * (ri-center)
            q = elements[abs(qi)]
            if qi < 0:
                q = f'"{q}"'
            print(f'  <text x="{ri[0]}" y="{ri[1]}" class="atnum">{q}</text>')

    if labels:
        print()
        for i, label in labels:
            ri = R[i]
            ri = afactor * (ri-center)
            print(f'  <text x="{ri[0]}" y="{ri[1]}" class="atnum">{label}</text>')


    if values:
        z_order = np.argsort([x[0][2] for x in values])
        maxval = max(map(lambda x: abs(x[-1]), values))
        for I in z_order:
            ri, value = values[I]
            value = value/maxval
            ri = afactor * (ri-center)
            mycolor = linear_grad(par.val_grad.gcolors, value)
            print(f'  <use fill="#{mycolor}" x="{ri[0]}" y="{ri[1]}" xlink:href="#values"/>')



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

    values = []
    for value in filter(lambda x: x[0]=='value', data):
        values.append(([*map(float, value[1:4])], float(value[4])))

    return np.array(Q), np.array(R), bonds, labels, values


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
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0x090909, 0xfcab28, 0xF0F000, 0x00F000, 0xA0A0A0,
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
      0xA0A0A0, 0xA0A0A0, 0xA0A0A0, 0x090909, 0xfcab28, 0xFFFF22, 0x22FF22, 0xA0A0A0,
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
      0x303030, 0x303030, 0x303030, 0x090909, 0x303030, 0x888800, 0x008800, 0x303030,
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
    elements = ['',
                'H', 'He',
                'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
                'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                ] + [str(i) for i in range(103, 128)]

    return np.array(radii), np.array(colors), np.array(colors_ini), np.array(colors_fin), np.array(elements)


if __name__=='__main__':
    radii, colors, colors_ini, colors_fin, elements = atom_parameters()
    parameters = parse_arguments(radii)
    if parameters.light_hydrogen:
         colors[1] = 0xd1d1d1
    Q, R, bonds, labels, values = mol_input()
    print_svg(Q, R, bonds, labels, values, radii, colors, colors_ini, colors_fin, elements, parameters)

