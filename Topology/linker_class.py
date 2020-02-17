import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class make_linked_proteins(object):


    def __init__(self, num_l1, num_l2, num_p, side):

        self.num_p = num_p
        self.num_l1 = num_l1
        self.num_l2 = num_l2



        self.sigma_proteins = 2.0
        self.sigma_linkers = 0.5

        self.side = side
        self.lattice_constant = 2.5*self.sigma_proteins
        self.Lx = self.side*self.lattice_constant
        self.Ly = self.side*self.lattice_constant
        self.Lz = 0.1

        self.lj_factor = 2**(1/6)
     

        self.coords = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles = ["\nAngles \n \n"]

        # m is number of atoms and k is number of 3 atom-molecules
        self.numAll = 0
        self.k = 0

        self.numBonds = 0
        self.bondId = 0

        self.numAngles = 0
        self.angleId = 0

        self.numTypes = 0
        self.type_mass_list =[]


        def init_lattce_2d(lat_con):
            num_x = int(self.Lx/lat_con)
            num_y = int(self.Ly/lat_con)

            lattice_out = np.zeros(shape=(num_x*num_y, 3))
            counter = 0
            for i in range(num_x):
                for j in range(num_y):

                    lattice_out[counter, 0] = lat_con * i - lat_con*num_x * 0.5
                    lattice_out[counter, 1] = lat_con * j - lat_con*num_y * 0.5
                    lattice_out[counter, 2] = 0.0
                    counter += 1

            return lattice_out



        self.lattice_sites = init_lattce_2d(self.lattice_constant)

        assert self.num_p + self.num_l1 + self.num_l2 <= self.lattice_sites.shape[0], "ERROR: not enough lattice sites, increase box size."



        lattice_inds = np.arange(self.lattice_sites.shape[0])
        self.start_inds_p = np.random.choice(lattice_inds, size = self.num_p, replace = False)

        self.lattice_inds_remain = np.array([item for item in lattice_inds if item not in set(list(self.start_inds_p))])


        self.start_inds_l1 = np.random.choice(self.lattice_inds_remain, size = self.num_l1, replace = False)

        self.lattice_inds_remain = np.array([item for item in lattice_inds if item not in set(list(self.start_inds_p) +  list(self.start_inds_l1))])

        self.start_inds_l2 = np.random.choice(self.lattice_inds_remain, size=self.num_l2, replace=False)


        # plt.scatter( self.lattice_sites[self.start_inds_l1][:, 0 ], 1+self.lattice_sites[self.start_inds_l1][:, 1 ])
        # plt.scatter( self.lattice_sites[self.start_inds_l2][:, 0 ], self.lattice_sites[self.start_inds_l2][:, 1 ], c ='g')
        # plt.scatter( self.lattice_sites[self.start_inds_p][:, 0 ], -1+self.lattice_sites[self.start_inds_p][:, 1 ], c ='r')
        # plt.show()
        #
        # import sys
        # sys.exit()
        test_list = list(self.start_inds_p) + list(self.start_inds_l1) + list(self.start_inds_l2)
        assert len(test_list) == self.num_l1 + self.num_l2 + self.num_p

    def make_LAT(self):

        particles_list = [1,2,3,4,5]
        particle_types = len(particles_list)
        self.numTypes += len(np.unique(particles_list))

        # for i in range(particle_types):
        #     to_add = [particles_list[i], 1]
        #     if to_add not in self.type_mass_list:
        #         self.type_mass_list.append(to_add)
        self.type_mass_list.append([1,1])
        self.type_mass_list.append([2,1])
        self.type_mass_list.append([3,1])
        self.type_mass_list.append([4,1])
        self.type_mass_list.append([5,1])

        for i in range(len(self.start_inds_p)):  ## go through chains
            indBuf = self.start_inds_p[i]

            self.k += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            # print(x,y,z)
            bonded_num = 0
            for n in range(5):
                self.numAll += 1
                if (n == 0):

                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " "+str(particles_list[n])+" 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")
                    bonded_num = self.numAll

                radial_dist = 1.0*self.lj_factor

                position_on_cirlce = np.array([[radial_dist,0],[0,-radial_dist],[-radial_dist,0], [0,radial_dist] ])
                if n > 0:

                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " +  str(particles_list[n]) + " 0 " + str(x - position_on_cirlce[n-1,0]) + " " + str(y - position_on_cirlce[n-1,1]) + " " + str(
                            z) + " 0 0 0 \n")

                    # self.bondId += 1
                    # self.bonds.append(
                    #     "\t " + str(self.bondId) + " 1 " + str(self.numAll) + " " + str(bonded_num) + " \n")
                    # self.numBonds += 1

                    # self.angleId += 1
                    # self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll - 1) + " " + str(
                    #     self.numAll - 2) + " " + str(self.numAll) + "\n")
                    # self.numAngles += 1



    def make_pLCgamma1(self):

        particles_list = [6, 7, 8, 9]  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += 4

        # self.type_mass_list = list(np.unique(self.type_mass_list))
        # for i in range(particle_types):
        #     self.type_mass_list.append([particles_list[i], 1.0])
        # self.type_mass_list = list(np.unique(self.type_mass_list))
        self.type_mass_list.append([6, 1])
        self.type_mass_list.append([7, 1])
        self.type_mass_list.append([8, 1])
        self.type_mass_list.append([9, 1])

        for i in range(len(self.start_inds_l1)):  ## go through chains
            indBuf = self.start_inds_l1[i]

            self.k += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            print(x,y,z)
            bonded_num = 0
            radial_dist = 1.0 * self.lj_factor
            position_on_cirlce = radial_dist*np.ones(shape = (3,2 ))
            ligand_angles = 2.0*np.pi*np.array([1/3, 2/3, 3/3])
            for i in range(len(ligand_angles)):
                position_on_cirlce[i, 0] = np.sin(ligand_angles[i])
                position_on_cirlce[i, 1] = np.cos(ligand_angles[i])


            for n in range(particle_types):
                self.numAll += 1
                if (n == 0):
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")
                    bonded_num = self.numAll



                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(
                            x - position_on_cirlce[n - 1, 0]) + " " + str(
                            y - position_on_cirlce[n - 1, 1]) + " " + str(
                            z) + " 0 0 0 \n")

                    # self.bondId += 1
                    # self.bonds.append(
                    #     "\t " + str(self.bondId) + " 1 " + str(self.numAll) + " " + str(
                    #         bonded_num) + " \n")
                    # self.numBonds += 1
                    #
                    # self.angleId += 1
                    # self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll - 1) + " " + str(
                    #     self.numAll - 2) + " " + str(self.numAll) + "\n")
                    # self.numAngles += 1

    def make_SOS1(self):

        particles_list = [10, 11, 12 ,13,14]  ## first one is central body
        particle_types = len(particles_list)
        self.numTypes += 5

        # self.type_mass_list = list(np.unique(self.type_mass_list))
        # for i in range(particle_types):
        #     self.type_mass_list.append([particles_list[i], 1.0])
        # self.type_mass_list = list(np.unique(self.type_mass_list))
        # self.type_mass_list.append([5, 1])
        self.type_mass_list.append([10, 1])
        self.type_mass_list.append([11, 1])
        self.type_mass_list.append([12, 1])
        self.type_mass_list.append([13, 1])
        self.type_mass_list.append([14, 1])

        for i in range(len(self.start_inds_l2)):  ## go through chains
            indBuf = self.start_inds_l2[i]

            self.k += 1
            x = self.lattice_sites[indBuf, 0]
            y = self.lattice_sites[indBuf, 1]
            z = self.lattice_sites[indBuf, 2]

            # print(x,y,z)
            bonded_num = 0
            for n in range(particle_types):
                self.numAll += 1
                if (n == 0):
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(x) + " " + str(y) + " " + str(
                            z) + " 0 0 0 \n")
                    bonded_num = self.numAll

                radial_dist = 1.0 * self.lj_factor
                position_on_cirlce = np.array(
                    [[radial_dist, 0], [-radial_dist, 0], [0, radial_dist], [0, -radial_dist]])
                if n > 0:
                    self.coords.append(
                        "\t " + str(self.numAll) + " " + str(self.k) + " " + str(
                            particles_list[n]) + " 0 " + str(
                            x - position_on_cirlce[n - 1, 0]) + " " + str(
                            y - position_on_cirlce[n - 1, 1]) + " " + str(
                            z) + " 0 0 0 \n")

                    # self.bondId += 1
                    # self.bonds.append(
                    #     "\t " + str(self.bondId) + " 1 " + str(self.numAll) + " " + str(
                    #         bonded_num) + " \n")
                    # self.numBonds += 1

                    # self.angleId += 1
                    # self.angles.append("\t " + str(self.angleId) + " 1 " + str(self.numAll - 1) + " " + str(
                    #     self.numAll - 2) + " " + str(self.numAll) + "\n")
                    # self.numAngles += 1



