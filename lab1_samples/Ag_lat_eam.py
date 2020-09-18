from labutil.src.plugins.lammps import *
from ase.spacegroup import crystal
from ase.build import *
import numpy
import matplotlib.pyplot as plt


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

# pair_style lj/cut 4.5
# pair_coeff 1 1 0.3450 2.644 4.5

# ---------- 3. Run single point calculation  ---------------------
thermo_style custom step pe lx ly lz press pxx pyy pzz
run 0

#-- include optimization of the unit cell parameter
#fix 1 all box/relax iso 0.0 vmax 0.001

#-- enable optimization of atomic positions (and the cell)
#min_style cg
#minimize 1e-10 1e-10 1000 10000


# ---- 4. Define and print useful variables -------------
variable natoms equal "count(all)"
variable totenergy equal "pe"
variable length equal "lx"

print "Total energy (eV) = ${totenergy}"
print "Number of atoms = ${natoms}"
print "Lattice constant (Angstoms) = ${length}"
        """

def make_struc(alat, supercell_size):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    multiplier = numpy.identity(3) * supercell_size
    ase_supercell = make_supercell(unitcell, multiplier)
    structure = Struc(ase2struc(ase_supercell))
    return structure

def make_struc_vac(alat, supercell_size):
    """
    Creates the fcc crystal structure with central atom missing using ASE.
    :param alat: Lattice parameter in angstrom
    :param supercell_size: supercell multiplier
    :return: structure object converted from ase
    """
    unitcell = crystal('Ag', [(0, 0, 0)], spacegroup=225, cellpar=[alat, alat, alat, 90, 90, 90])
    multiplier = numpy.identity(3) * supercell_size
    ase_supercell = make_supercell(unitcell, multiplier)
    for i in range(len(ase_supercell.positions)):
        if (ase_supercell.positions[i] == [supercell_size*alat/2, supercell_size*alat/2, supercell_size*alat/2]).all():
            ase_supercell.pop(i)
            break
    structure = Struc(ase2struc(ase_supercell))
    return structure


def compute_energy(alat, template):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Ag_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Ag"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))
    struc = make_struc(alat=alat)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)
    return energy, lattice

def compute_vac_energy(alat, supercell_size, template):
    """
    Make an input template and select potential and structure, and the path where to run
    Computes energy of supercell with and without vacancy
    """
    potpath = os.path.join(os.environ['LAMMPS_POTENTIALS'],'Ag_u3.eam')
    potential = ClassicalPotential(path=potpath, ptype='eam', element=["Ag"])
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab1", str(alat)))

    struc = make_struc(alat=alat, supercell_size=supercell_size)
    output_file = lammps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy, lattice = get_lammps_energy(outfile=output_file)

    struc_vac = make_struc_vac(alat=alat, supercell_size=supercell_size)
    output_file = lampmps_run(struc=struc, runpath=runpath, potential=potential, intemplate=template, inparam={})
    energy_vac, lattice_vac = get_lammps_energy(outfile=output_file)

    print(f'Energy: {energy}')
    print(f'Energy with vacancy: {energy_vac}')
    return energy, energy_vac


def lattice_scan():
    alat_list = numpy.linspace(3.8, 4.3, 6)
    energy_list = [compute_energy(alat=a, template=input_template)[0] for a in alat_list]
    print(alat_list)
    print(energy_list)
    #plt.plot(alat_list, energy_list)
    #plt.show()


if __name__ == '__main__':
    # put here the function that you actually want to run
    #lattice_scan()
    #compute_energy(4.1, input_template)
    compute_vac_energy(4.08999994375577, input_template)
