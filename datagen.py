import numpy as np
import math
import random

Numatoms, Numbonds, Numangles, Numdihedrals, Nimpropers = 0, 0, 0, 0, 0
Numatomtype, Numbondtype, Numangletype, Numdihedraltype, Numpropertype = 5,1,1,0,0
BoxX, BoxY, BoxZ = 100, 100, 100
Mass = 1.0

#this is random decision function
def decision(probability):
    return random.random() < probability

#write the head of data file
def write_head(Numatoms, Numbonds, Numangles, Numdihedrals, Nimpropers,Numatomtype, \
Numbondtype, Numangletype, Numdihedraltype, Numpropertype,BoxX, BoxY, BoxZ):
    out_string ='''LAMMPS Membrane System

       {0}  atoms
       {1}  bonds
       {2}  angles
       {3}  dihedrals
       {4}  impropers

       {5}  atom types
       {6}  bond types
       {7}  angle types
       {8}  dihedral types
       {9}  improper types\n'''\
       .format(Numatoms, Numbonds, Numangles, Numdihedrals, Nimpropers,\
       Numatomtype, Numbondtype, Numangletype, Numdihedraltype, Numpropertype)
    out_string +='''
        0.000000  {0} xlo xhi
        0.000000  {1} ylo yhi
        0.000000  {2} zlo zhi\n'''.format("%.6f" %BoxX, "%.6f" %BoxY, "%.6f" %BoxZ)
    f.write(out_string)

#write the masses information to data file
def write_mass(Mass):

    out_string = "\nMasses\n"+'\n'
    f.write(out_string)
    for i in range(Numatomtype):
        out_string  = '{0:{fill}{align}6}'.format( i+1 , fill=' ', align='>')
        out_string += '{0:{fill}{align}6}'.format( Mass , fill=' ', align='>')
        out_string += '\n'
        print(out_string)
        f.write(out_string)

#write the information of atoms to the data file
# use shape function to generate specific coordinates
def write_atoms(BoxX,BoxY,BoxZ):
    out_string = "\nAtoms\n"
    out_string +="  # atom_id, molecule_id, atom_type, charge, x, y, z \n"
    f.write(out_string)
#    membrane(BoxX,BoxY,BoxZ)
    return_id = membrane(BoxX,BoxY,BoxZ)
#    return_id = [0,0]
    atom_id = return_id[0]
    atom_id_lipid = atom_id
    molecule_id = return_id[1]
    atom_id = sphere(atom_id,molecule_id)
    return (atom_id, atom_id_lipid)



def membrane(BoxX,BoxY,BoxZ):
    Sigma = 1
    Nx = int(BoxX/Sigma)
    Ny = int(BoxY//Sigma)
    Nz = int(6)
    atom_id = 0
    molecule_id = 0
    for i in range(Nx):
        y = 0.5 * Sigma
        x = Sigma*i+0.5 * Sigma
        for j in range(Ny):
            z = BoxZ/2 + 3 * Sigma
            y = Sigma*j+0.5 * Sigma
            molecule_id += 1
            for k in range(Nz):
                if k<=2:
                    z = BoxZ/2 + 3 * Sigma-Sigma*k
                else:
                    molecule_id += 1-math.ceil((k-3)/2)
		    z = BoxZ/2 - 3 * Sigma+Sigma*(k-2)
                atom_id +=1
                if k%3==0:
                    if decision(0.5):
                        atom_type = 2
                    else:
                        atom_type = 1
                else:
                    atom_type = 3
		molecule_id = int(molecule_id)
                writedata(atom_id,molecule_id,atom_type,0,x,y,z)
    return (atom_id,molecule_id)

def writedata(atom_id,molecule_id,atom_type,charge,x,y,z):
    out_string = '{0:{fill}{align}10}'.format(atom_id , fill=' ', align='>')
    out_string += '{0:{fill}{align}10}'.format(molecule_id , fill=' ', align='>')
    out_string += '{0:{fill}{align}10}'.format(atom_type , fill=' ', align='>')
    out_string += '{0:{fill}{align}12}'.format("%.6f" % charge , fill=' ', align='>')
    out_string += '{0:{fill}{align}12}'.format("%.6f" % x , fill=' ', align='>')
    out_string += '{0:{fill}{align}12}'.format("%.6f" % y , fill=' ', align='>')
    out_string += '{0:{fill}{align}12}'.format("%.6f" % z , fill=' ', align='>')
    out_string += '\n'
    f.write(out_string)

def sphere(atom_id,molecule_id):
    Sigma = 1
    R = 5
    center = [50,50,80]
    atom_id = atom_id
    molecule_id = molecule_id+1
    top_id =atom_id+1
    bottom_id =atom_id+2
    writedata(top_id,molecule_id,4,0,center[0],center[1],center[2]+R)
    writedata(bottom_id,molecule_id,5,0,center[0],center[1],center[2]-R)
    atom_id = atom_id+2
    dtheta = 2*np.arctan(0.5*Sigma/R)
    ntheta = int(math.ceil(np.pi/dtheta)-2)
    for i in range(1,ntheta+1):
        theta = np.pi*i/ntheta
        r = R*np.sin(theta)
        dphi  = 2*np.arctan(0.5*Sigma/r)
        nphi = int(math.ceil(np.pi/dphi)*2)
        for j in range(1,nphi+1):
            phi  = 2*np.pi*j/nphi
            x=R*np.sin(theta)*np.cos(phi)+center[0]
            y=R*np.sin(theta)*np.sin(phi)+center[1]
            z=R*np.cos(theta)+center[2]
            atom_id += 1
            if decision(0.5):
                atom_type = 4
            else:
                atom_type = 5
            writedata(atom_id,molecule_id,atom_type,0,x,y,z)
    return atom_id




def write_bond(Numatoms_lipid):
    out_string = "\nBonds\n"
    out_string +="  # bond_id, bond_type, atom1_id, atom2_id \n"
    f.write(out_string)
    bond_id = 0
    bond_type = 1
    for i in range(0,Numatoms_lipid,3):
        for j in range(2):
            atom1_id = i+j+1
            atom2_id = i+j+2
            bond_id +=1
            out_string = '{0:{fill}{align}8}'.format(bond_id , fill=' ', align='>')
            out_string += '{0:{fill}{align}8}'.format(bond_type , fill=' ', align='>')
            out_string += '{0:{fill}{align}8}'.format(atom1_id , fill=' ', align='>')
            out_string += '{0:{fill}{align}8}'.format(atom2_id , fill=' ', align='>')
            out_string += '\n'
            f.write(out_string)
    return bond_id

def write_angle(Numatoms_lipid):
    out_string = "\nAngles\n"
    out_string +="  # angle_id, angle_type, atom1_id, atom2_id , atom3_id\n"
    f.write(out_string)
    angle_id = 0
    angle_type = 1
    for i in range(0,Numatoms_lipid,3):
        atom1_id = i+1
        atom2_id = i+2
        atom3_id = i+3
        angle_id +=1
        out_string = '{0:{fill}{align}8}'.format(angle_id , fill=' ', align='>')
        out_string += '{0:{fill}{align}8}'.format(angle_type , fill=' ', align='>')
        out_string += '{0:{fill}{align}8}'.format(atom1_id , fill=' ', align='>')
        out_string += '{0:{fill}{align}8}'.format(atom2_id , fill=' ', align='>')
        out_string += '{0:{fill}{align}8}'.format(atom3_id , fill=' ', align='>')
        out_string += '\n'
        f.write(out_string)
    return angle_id





f = open("temp.delete" ,"w+")
Numatoms_array = write_atoms(BoxX,BoxY,BoxZ)
Numatoms = Numatoms_array[0]
Numatoms_lipid = Numatoms_array[1]
Numbonds = write_bond(Numatoms_lipid)
Numangles = write_angle(Numatoms_lipid)
f.close()

f = open("data.system" ,"w+")
print(Numatoms)
write_head(Numatoms, Numbonds, Numangles, Numdihedrals, Nimpropers,Numatomtype, \
Numbondtype, Numangletype, Numdihedraltype, Numpropertype,BoxX, BoxY, BoxZ)
write_mass(Mass)
write_atoms(BoxX,BoxY,BoxZ)
write_bond(Numatoms_lipid)
write_angle(Numatoms_lipid)
f.close()
