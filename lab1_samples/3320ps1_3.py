#from labutil.src.plugins.lammps import *
from ase.build import *
from ase.io import write

input_template = """
# ---------- 1. Initialize simulation ---------------------
units metal
atom_style atomic
dimension  3
boundary   p p p
read_data $DATAINPUT

# ---------- 2. Specify interatomic potential ---------------------
pair_style eam
pair_coeff * * $POTENTIAL

#pair_style lj/cut 4.5
#pair_coeff 1 1 0.3450 2.644 4.5

# ----------pre 3 Dump attempt ---------------
dump myDump all atom 100 dump.atom

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
#run 0

#-- include optimization of the unit cell parameter
fix 1 all box/relax iso 0.0 vmax 0.001

#-- enable optimization of atomic positions (and the cell)
min_style cg
minimize 1e-10 1e-10 1000 10000


# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
        """

def make_slab(alat, size_tuple):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    slab = fcc100('Ag', a=alat, size=size_tuple, vacuum=10.0)
    write('slab.cif', slab)
    #structure = Struc(ase2struc(slab))
    structure =3
    return structure

def compute_energy(alat, size_tuple, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Ag_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Ag"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))
    struc = make_slab(alat, size_tuple)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)

    print(f'Cohesive Energy: {energy}')

    return energy, lattice

if __name__ == '__main__':
    # put here the function that you actually want to run
    make_slab(4.089999944, (4,4,20))
    #compute_energy(4.089999944, (4, 4, 20), input_template)
