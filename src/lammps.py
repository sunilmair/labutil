from string import Template
from labutil.src.util import *
from labutil.src.objects import *
import os
from ase.io.lammpsrun import read_lammps_dump


def write_lammps_data(struc, runpath):
    """Make LAMMPS struc data"""
    datatxt = 'Header of the LAMMPS data file \n\n'
    datatxt += '{} atoms \n'.format(struc.n_atoms)
    datatxt += '{} atom types\n\n'.format(struc.n_species)

    datatxt += '0.0 {}  xlo xhi\n'.format(struc.cell[0][0])
    datatxt += '0.0 {}  ylo yhi\n'.format(struc.cell[1][1])
    datatxt += '0.0 {}  zlo zhi\n'.format(struc.cell[2][2])

    datatxt += '\nMasses \n\n'
    for sym, sp in struc.species.items():
        datatxt += '{}  {}\n'.format(sp['kind'], sp['mass'])
    # Write atom positions in angstrom
    datatxt += '\nAtoms # atomic \n\n'
    for index, site in enumerate(struc.sites):
        datatxt += '{} {} {:1.5f} {:1.5f} {:1.5f} \n'.format(index + 1, struc.species[site[0]]['kind'], *site[1])

    datafile = os.path.join(runpath.path, 'lammps.data')
    write_file(datafile, datatxt)
    return File({'path': datafile})


def write_lammps_input(datafile, runpath, in_template, potential=None):
    """make Lammps input script"""
    if potential:
        ppath = potential.path
    else:
        ppath = ''

    dumpfile = os.path.join(runpath.path, 'lammps.dump')
    subst = {'DATAINPUT': datafile.path, 'POTENTIAL': ppath, 'DUMP': dumpfile}
    inptxt = Template(in_template).safe_substitute(subst)
    infile = os.path.join(runpath.path, 'lammps.in')
    write_file(infile, inptxt)
    return File({'path': infile})


def lammps_run(struc, runpath, in_template, potential=None):
    """
    Run the LAMMPS code with provided inputs
    """
    lammps_code = ExternalCode({'path': os.environ['LAMMPS_COMMAND']})
    prepare_dir(runpath.path)
    datafile = write_lammps_data(struc=struc, runpath=runpath)
    infile = write_lammps_input(datafile=datafile, potential=potential, runpath=runpath, in_template=in_template)
    logfile = File({'path': os.path.join(runpath.path, 'lammps.log')})
    outfile = File({'path': os.path.join(runpath.path, 'lammps.out')})
    lammps_command = "{} -in {} -log {} > {}".format(lammps_code.path, infile.path, logfile.path, outfile.path)
    run_command(lammps_command)
    return outfile


def get_lammps_energy(outfile):
    """
    Read the final lattice parameter and energy from output file
    """
    energy = None
    lattice = None
    with open(outfile.path, 'r') as fout:
        for line in fout.readlines():
            if 'Total energy' in line:
                energy = float(line.split(" = ")[1])
            if 'Lattice constant (Angstoms)' in line:
                lattice = float(line.split(" = ")[1])
    return energy, lattice


def parse_structure_dump(runpath, dumpfilename):
    """
    Using ASE to parse the LAMMPS output. Needs testing.
    """
    dumpfile = os.path.join(runpath.path, dumpfilename)
    last_structure = read_lammps_dump(dumpfile)
    strucfile = os.path.join(runpath.path, 'struc.cif')
    ase.io.write(strucfile, last_structure)
