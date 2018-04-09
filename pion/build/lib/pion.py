#######################################################################################################################
# pion@Hadron is open-sourced program based on python. So, you can freely adjust or fix this code. If you have some   #
#problem or if you have done awesome modification to this program, please contact me: Hadron@yonsei.ac.kr . I appreci-#
#ate to everyone to help me.                                                                                          #
#                    Atom color scheme is used at VMD's(visualized molecular dynamics) one.                           #
#######################################################################################################################

from mayavi import mlab
import numpy as np
import argparse
import atom_data as ad                              #Import atom database: atom_data

#####################################
#         Argument Module           #
#####################################

parser = argparse.ArgumentParser(description="pion@Hadron cube file visualization program",
                                 epilog="Hadron@Yonsei TCCL 2017")
parser.add_argument('-i', help="Input Cube File")
# parser.add_argument('-o', help='image file output')
parser.add_argument('-b', help="Bond Information File, Automatically Ignored when -a False", default="bond.dat")
parser.add_argument('-bt', help="Bond Information File Format (ascii or mol2), default=ascii", default="ascii")
parser.add_argument('-v', help="Show version information", action='version',
                    version='pion@Hadron alpha 1.00')
parser.add_argument('-a', help="Only display Atoms without Bond, default=True", type=bool, default=True)
parser.add_argument('-l', help="minimum transparent of density", type=float, default = 0.6)
parser.add_argument('-t', help='Maximum transparent of density', type=float, default = 1.0)
parser.add_argument('-r', help='Set Atom Resolution', type=int, default = 100)


args = parser.parse_args()

print("      3.141592653589793238462643383279  ")
print("    5028841971693993751058209749445923  ")
print("   07816406286208998628034825342117067  ")
print("   9821    48086         5132           ")
print("  823      06647        09384           ")
print(" 46        09550        58223           ")
print(" 17        25359        4081            ")
print("           2848         1117            ")
print("           4502         8410            ")
print("           2701         9385            ")
print("          21105        55964            ")
print("          46229        48954            ")
print("          9303         81964            ")
print("          4288         10975            ")
print("         66593         34461            ")
print("        284756         48233            ")
print("        78678          31652        71  ")
print("       2019091         456485       66  ")
print("      9234603           48610454326648  ")
print("     2133936            0726024914127   ")
print("     3724587             00660631558    ")
print("     817488               152092096     ")


#####################################
#       Read Data from input        #
#####################################

cube = open(args.i, 'rb')

for k in range(2):
    line = cube.readline()

(atom_number, origin_x, origin_y, origin_z) = cube.readline().split()
(x_grid, x_pixel, zero, zero2) = cube.readline().split()
(y_grid, zero, y_pixel, zero2) = cube.readline().split()
(z_grid, zero, zero2, z_pixel) = cube.readline().split()

mol_atoms = []
atoms_x = []
atoms_y = []
atoms_z = []

mlab.figure(1, bgcolor=(0, 0, 0), size=(400, 400))
mlab.clf()

for i in range(1, int(atom_number)+1):
    [number, charge, x_coord, y_coord, z_coord] = cube.readline().split()
    mol_atoms.append(min(int(number), 55))
    atoms_x.append(float(x_coord)-float(origin_x))
    atoms_y.append(float(y_coord)-float(origin_y))
    atoms_z.append(float(z_coord)-float(origin_z))

atoms_x = np.array(atoms_x) / float(x_pixel) + 1.0
atoms_y = np.array(atoms_y) / float(y_pixel) + 1.0
atoms_z = np.array(atoms_z)/ float(z_pixel) + 1.0

for i in range(0, int(atom_number)):
    mlab.points3d(atoms_x[i], atoms_y[i], atoms_z[i],
                  scale_factor=ad.atom_scale[mol_atoms[i]],
                  color=ad.atom_color[mol_atoms[i]-1],
                  resolution=args.r,
                  scale_mode='none')

if args.a == True:
    bond = open(args.b, 'rb')
    while True:
        line = bond.readline()
        if not line:
            break
        [atom1, atom2] = line.split()
        mlab.plot3d([atoms_x[int(atom1) - 1], atoms_x[int(atom2) - 1]],
                    [atoms_y[int(atom1) - 1], atoms_y[int(atom2) - 1]],
                    [atoms_z[int(atom1) - 1], atoms_z[int(atom2) - 1]],
                    [1, 2], tube_radius=0.4, color=(1,1,1))

str = ' '.join(file(args.i).readlines()[6+int(atom_number):])
#data = np.fromstring(str, sep=' ').reshape((40, 40, 40))
data = np.fromstring(str, sep=' ').reshape((int(x_grid), int(y_grid), int(z_grid)))

source = mlab.pipeline.scalar_field(data/data.max())
min = data.min()
max = data.max()
vol = mlab.pipeline.volume(source, vmin=args.l, vmax=args.t)
#vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
#                           vmax=min + 0.9 * (max - min))

#mlab.axes()

mlab.view(132, 54, 45, [np.mean(atoms_x), np.mean(atoms_y), np.mean(atoms_z)])

print('       _               ____                  _                 ')
print(' _ __ (_) ___  _ __   / __ \  /\  /\__ _  __| |_ __ ___  _ __  ')  
print('| \'_ \| |/ _ \| \'_ \ / / _` |/ /_/ / _` |/ _` | \'__/ _ \| \'_ \ ') 
print('| |_) | | (_) | | | | | (_| / __  / (_| | (_| | | | (_) | | | |') 
print('| .__/|_|\___/|_| |_|\ \__,_\/ /_/ \__,_|\__,_|_|  \___/|_| |_|')  
print('|_|                   \____/                                   ') 


#print("            **                  ")
#print("    ******                      ")
#print("    **   ** **  ******  ******* ")
#print("    **   ** ** **    **  **   **")
#print("    ******  ** **    **  **   **")
#print("    **      ** **    **  **   **")
#print("    **      **  ******  ***   **")

mlab.show()
