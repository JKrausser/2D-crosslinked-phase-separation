import argparse
from Topology.linker_class import make_linked_proteins
import itertools
import os
import numpy as np

def write_in_script( eps_pp, eps_pl, eps_ll, real,folderPattern,filePattern):

    filename = "Input/Scripts/Input_%s_%i.in" % (filePattern, real)

    f = open(filename, "w")

    lj_factor = 2**(1/6)

    global_cutoff = 3.0
    protein_sigma = 2.0
    linker_sigma = 0.5
    cross_sigma = 0.5 * (protein_sigma + linker_sigma)



    f.write("units                  lj                                      \n")
    f.write("dimension 2\n")
    f.write("atom_style             full  \n")
    f.write("boundary               p p p   \n")
    f.write("read_data              %s \n"%configName)




    f.write("neighbor               0.3 bin\n")
    f.write("neigh_modify           every 1 delay 1\n")

    kBond_protein = 100.0
    kBond_linker = 10.0
    kBend = 10
    f.write("bond_style             harmonic\n")
    f.write(f"bond_coeff             1 {kBond_protein} {0.5*protein_sigma*lj_factor}\n")
    f.write(f"bond_coeff             2 {kBond_linker} {linker_sigma*lj_factor}\n")
    f.write("angle_style            harmonic\n")
    # f.write(f"angle_coeff            1 {kBend} 180\n")
    # f.write("angle_coeff            2 %lf 180\n"%kBend)


    f.write(
    '''
    ##======================================##
    ## Interactions
    ##======================================##
    \n''')

    eps_pp = 1
    eps_pl = 1
    eps_ll = 2
    ligand_range = 2
    f.write("pair_style             cosine/squared %lf \n"%global_cutoff)

    f.write("pair_coeff             * * %lf %lf %lf wca \n" % (0, lj_factor , lj_factor))  ## repulsive only



    repulsive_only = [1, 6, 10]
    for id in repulsive_only:
        for id2 in repulsive_only:
            if id2 >= id:
                f.write(f"pair_coeff             {id} {id2} %lf %lf %lf wca \n" % (eps_pp, lj_factor * protein_sigma, lj_factor * protein_sigma))  ## repulsive only

    # weak ligand interaction
    # lat_int_set = [2, 3, 4, 5]
    # plcgamma1_int_set = [7, 8, 9]
    # for id_lat in lat_int_set:
    #     for id_plc in plcgamma1_int_set:
    #         f.write(f"pair_coeff             {id_lat} {id_plc} {eps_ll} {lj_factor*linker_sigma} {ligand_range + lj_factor*linker_sigma} wca \n") ## repulsive only

    f.write(f"pair_coeff             2 9 {3*eps_ll} {lj_factor * linker_sigma} {ligand_range + lj_factor * linker_sigma} wca \n")  ## repulsive only
    f.write(f"pair_coeff             3 7 {3*eps_ll} {lj_factor * linker_sigma} {ligand_range + lj_factor * linker_sigma} wca \n")  ## repulsive only


    #SOS1 interaction
    sos_int_set = [11, 12, 13, 14]
    for id_sos in sos_int_set:
        f.write(f"pair_coeff             8 {id_sos} {eps_ll} {lj_factor * linker_sigma} {ligand_range + lj_factor * linker_sigma} wca \n")  ## repulsive only

    f.write(
    '''
    ##======================================##
    ## Setting up groups
    ##======================================##
    \n''')
    # f.write("group                  head type 1\n")


    f.write(
    '''
    ##======================================##
    ## Computes
    ##======================================##
    \n''')






    f.write(
    '''
    ##======================================##
    ## Fixes
    ##======================================##
    \n''')

    f.write("velocity               all create 1.0 1 \n")


    f.write("log                    %s/Log_%s_%i.dat\n"%(folderPattern,filePattern,real))
    


    f.write("fix                    fLANG all langevin 1.0 1.0 1.0 %i\n"%seed)
    # f.write("fix                    fNVE all nve\n")
    if rigid_molecules:
        f.write("neigh_modify           exclude molecule/intra all\n")
        f.write("fix                    rigidNVE all rigid/nve molecule\n")







    f.write(
    '''
    ##======================================##
    ## Output
    ##======================================##
    \n''')

    f.write("thermo                 %i\n"%dump_step)
    f.write("thermo_style           custom step temp press etotal epair\n")

    f.write("thermo_modify          flush yes\n")
    f.write("timestep               0.01\n")

    f.write(
    '''
    ##======================================##
    ## Equilibration before output
    ##======================================##
    \n''')



    f.write(
    '''
    ##======================================##
    ## Run with output
    ##======================================##
    \n''')


    f.write("fix enforce_2d all enforce2d\n")
    f.write("dump                       1 all custom %i %s/Movie_%s_%i.xyz id type mol x y z \n"%(dump_step, folderPattern, filePattern,real))

    f.write("run                    %i\n"%run_steps)
    f.close()





    




if __name__ == "__main__":

    rigid_molecules = True

    sim_seed = 5
    dump_step = 1e3
    eq_steps = 1e3
    run_steps = 1e6

    ##======================================##
    ## Simulation parameters
    ##======================================##

    eps_pp = 1.0 ## strength tip-tip interaction between proteins
    eps_pl = 1.0 ## strength tip-tip interaction between proteins
    eps_ll = 1.0 ## strength tip-tip interaction between proteins
    
    num_pLCgamma1 = 200 ##linker concentration
    num_SOS1 = 0 ##linker concentration
    num_LAT = 100 ## protein concentration
    side_box = 20

    parser = argparse.ArgumentParser(description="Script for generating membrane.")
    parser.add_argument('--eps_pp', '-eps_pp', dest='eps_pp', action='store', type=float, default=eps_pp, help='range of attraction')
    parser.add_argument('--eps_pl', '-eps_pl', dest='eps_pl', action='store', type=float, default=eps_pl, help='range of attraction')
    parser.add_argument('--eps_ll', '-eps_ll', dest='eps_ll', action='store', type=float, default=eps_ll, help='range of attraction')

    parser.add_argument('--num_l1', '-num_l1', dest='num_l1', action='store', type=int, default=num_pLCgamma1, help='range of attraction')
    parser.add_argument('--num_l2', '-num_l2', dest='num_l2', action='store', type=int, default=num_SOS1, help='range of attraction')
    parser.add_argument('--num_p', '-num_p', dest='num_p', action='store', type=int, default=num_LAT, help='range of attraction')


    parser.add_argument('-real','--real', dest='real', action='store', type=int, default=0, help='number of realisation (statistics)')
    parser.add_argument('-seed','--seed', dest='seed', action='store', type=int, default=100, help='number of realisation (statistics)')

    args = parser.parse_args()

    eps_pp = args.eps_pp
    eps_pl = args.eps_pl
    eps_ll = args.eps_ll

    num_l1 = args.num_l1
    num_l2 = args.num_l2
    num_p = args.num_p
    


    real = args.real
    seed = args.seed


    
    system = make_linked_proteins(num_l1, num_l2, num_p,side_box)
    system.make_LAT()
    system.make_pLCgamma1()
    system.make_SOS1()


    filePattern = f"2d_linker_eps_pp{eps_pp}_eps_pl{eps_pl}_eps_ll{eps_ll}_n_l1{num_l1}_n_l2{num_l2}_n_p{num_p}_{real}"
    folderPattern = f"Results_n_l1{num_l1}_n_l2{num_l2}_n_p{num_p}"

    if not os.path.exists(folderPattern):
        os.makedirs(folderPattern)

    configName = "Input/Configuration/Config_%s.dat" % filePattern

    header = ["LAMMPS Description \n \n",
              "\t " + str(system.numAll) + " atoms \n \t " + str(system.numBonds) +
              " bonds \n \t " + str(system.numAngles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t "+str(system.numTypes)+" atom types \n \t 2 bond types \n \t 0 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(-system.Lx*0.5) + " " + str(system.Lx*0.5) + " xlo xhi\n \t", str(-system.Ly*0.5) + " " + str(system.Ly*0.5) + " ylo yhi \n \t",
              str(-system.Lz*0.5) + " " + str(system.Lz*0.5) + " zlo zhi\n"]

    header.append("\nMasses \n \n")
    for i in range(len(system.type_mass_list)):
        header.append("\t %i %.4f \n"%(system.type_mass_list[i][0],system.type_mass_list[i][1]))

    f = open(configName, "w")

    for item in header:
        f.write("%s " % item)

    for item in system.coords:
        f.write("%s " % item)

    if not rigid_molecules:
        for item in system.bonds:
            f.write("%s" % item)

    # for item in system.angles:
    #     f.write("%s" % item)

    f.close()
    write_in_script(eps_pp, eps_pl, eps_ll, real, folderPattern,filePattern)
















