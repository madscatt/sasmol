'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#	SASPDBRX
#
#	1/26/2012	--	initial coding				: hz/jc
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
import sasmol as sasmol

"""
The topology dictionary looks like the following:

pprint.pprint(o.topology_info['NTER'],width=100)
{'ATOM': [['N', 'NH3', '-0.30'],
          ['HT1', 'HC', '0.33'],
          ['HT2', 'HC', '0.33'],
          ['HT3', 'HC', '0.33'],
          ['CA', 'CT1', '0.21'],
          ['HA', 'HB', '0.10']],
 'BOND': [['HT1', 'N'], ['HT2', 'N'], ['HT3', 'N']],
 'DELE': [['ATOM', 'HN']],
 'DONO': ['HT1', 'N', 'HT2', 'N', 'HT3', 'N'],
 'IC': [['HT1', 'N', 'CA', 'C', '0.0000', '0.0000', '180.0000', '0.0000', '0.0000'],
        ['HT2', 'CA', '*N', 'HT1', '0.0000',
            '0.0000', '120.0000', '0.0000', '0.0000'],
        ['HT3', 'CA', '*N', 'HT2', '0.0000', '0.0000', '120.0000', '0.0000', '0.0000']],
 'TOTAL_CHARGE': '1.00'}

pprint.pprint(o.topology_info['ALA'],width=100)
{'ACCE': ['O', 'C'],
 'ATOM': [['N', 'NH1', '-0.47'],
          ['HN', 'H', '0.31'],
          ['CA', 'CT1', '0.07'],
          ['HA', 'HB', '0.09'],
          ['CB', 'CT3', '-0.27'],
          ['HB1', 'HA', '0.09'],
          ['HB2', 'HA', '0.09'],
          ['HB3', 'HA', '0.09'],
          ['C', 'C', '0.51'],
          ['O', 'O', '-0.51']],
 'BOND': [['CB', 'CA'],
          ['N', 'HN'],
          ['N', 'CA'],
          ['C', 'CA'],
          ['C', '+N'],
          ['CA', 'HA'],
          ['CB', 'HB1'],
          ['CB', 'HB2'],
          ['CB', 'HB3']],
 'DONO': ['HN', 'N'],
 'DOUB': [['O', 'C']],
 'IC': [['-C', 'CA', '*N', 'HN', '1.3551', '126.4900', '180.0000', '115.4200', '0.9996'],
        ['-C', 'N', 'CA', 'C', '1.3551', '126.4900',
            '180.0000', '114.4400', '1.5390'],
        ['N', 'CA', 'C', '+N', '1.4592', '114.4400',
            '180.0000', '116.8400', '1.3558'],
        ['+N', 'CA', '*C', 'O', '1.3558', '116.8400',
            '180.0000', '122.5200', '1.2297'],
        ['CA', 'C', '+N', '+CA', '1.5390', '116.8400',
            '180.0000', '126.7700', '1.4613'],
        ['N', 'C', '*CA', 'CB', '1.4592', '114.4400',
            '123.2300', '111.0900', '1.5461'],
        ['N', 'C', '*CA', 'HA', '1.4592', '114.4400',
            '-120.4500', '106.3900', '1.0840'],
        ['C', 'CA', 'CB', 'HB1', '1.5390', '111.0900',
            '177.2500', '109.6000', '1.1109'],
        ['HB1', 'CA', '*CB', 'HB2', '1.1109', '109.6000',
            '119.1300', '111.0500', '1.1119'],
        ['HB1', 'CA', '*CB', 'HB3', '1.1109', '109.6000', '-119.5800', '111.6100', '1.1114']],
 'IMPR': [['N', '-C', 'CA', 'HN'], ['C', 'CA', '+N', 'O']],
"""

import os
import numpy
import copy
import sasconfig as sasconfig

class Topology(object):

    '''
    This class contains charmm topolgy information used other modules.
    '''
    
    def __init__(self):
        pass

    def add(self, dictionary, key, value):
        '''
        This method is only used by read_topolgy method
        It will check whether the key is in the dictionary:
        if yes, append value to dictionary[key]
        if no, initialize dictionary[key] as [value]
        '''
        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)

    def read_charmm_topology(self, topology_file_path='', topology_file_name='top_all27_prot_na.inp'):
        '''
        Read and parse the charmm topology file
        A comprehensive dictionary (topology_info) will be built to store all the topology information
        The present strategy for parsing topology file is by splitting the words of each line
        '''
        error = []
        self.topology_info = {}
        lines = open(os.path.join(topology_file_path,
                                  topology_file_name), 'r').readlines()
        for line in lines:
            line = line.strip()
            if (len(line) > 0 and line[0] != '!'):
                words = line.split()
                if words[0][0:4] == 'MASS':
                    self.add(self.topology_info, 'MASS', words[1:4])
                #
                elif words[0][0:4] == 'DECL':
                    for ind in range(1, 2):
                        self.add(self.topology_info, 'DECL', words[ind])
                #
                elif words[0][0:4] == 'DEFA':  # Need further parsing
                    for ind in range(1, 5):
                        self.add(self.topology_info, 'DEFA', words[ind])
                #
                elif words[0][0:4] == 'AUTO':  # Need further parsing
                    for ind in range(1, 3):
                        self.add(self.topology_info, 'AUTO', words[ind])
                #
                elif words[0][0:4] == 'RESI' or words[0][0:4] == 'PRES':
                    cur_res = words[1]
                    self.topology_info[cur_res] = {}
                    self.topology_info[cur_res]['TOTAL_CHARGE'] = words[2]
                #
                else:
                    # elif 'cur_res' in locals():
                    #
                    if words[0] == 'ATOM':
                        self.add(self.topology_info[
                                 cur_res], 'ATOM', words[1:4])
                    #
                    elif words[0][0:4] == 'BOND':
                        ind = 1
                        while (ind + 2 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the BOND line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'BOND', words[ind:ind + 2])
                            ind += 2
                    #
                    elif words[0][0:4] == 'DOUB':
                        ind = 1
                        while (ind + 2 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DOUBLE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DOUB', words[ind:ind + 2])
                            ind += 2
                    #
                    elif words[0][0:4] == 'IMPR':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the IMPR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'IMPR', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'CMAP':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the CMAP line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'CMAP', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'ANGL':
                        ind = 1
                        while (ind + 3 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the ANGLE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'ANGL', words[ind:ind + 3])
                            ind += 3
                    #
                    elif words[0][0:4] == 'THET':
                        ind = 1
                        while (ind + 3 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the THET line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'THET', words[ind:ind + 3])
                            ind += 3
                    #
                    elif words[0][0:4] == 'DIHE':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DIHE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DIHE', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'DONO':  # Need to check
                        ind = 1
                        while (ind + 1 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DONOR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DONO', words[ind])
                            ind += 1
                    #
                    elif words[0][0:4] == 'ACCE':  # Need to check
                        ind = 1
                        while (ind + 1 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the ACCEPTOR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'ACCE', words[ind])
                            ind += 1
                    #
                    elif words[0][0:4] == 'IC':
                        self.add(self.topology_info[
                                 cur_res], 'IC', words[1:10])
                    #
                    elif words[0][0:4] == 'BUIL':
                        self.add(self.topology_info[
                                 cur_res], 'BUIL', words[1:10])
                    #
                    elif words[0][0:4] == 'PATC':  # Need further parsing
                        self.add(self.topology_info[
                                 cur_res], 'PATC', words[1:])
                    #
                    elif words[0][0:4] == 'DELE':
                        if 'DELE' not in self.topology_info[cur_res]:
                            self.topology_info[cur_res]['DELE'] = {}
                        prop = words[1][0:4]
                        if prop == 'ATOM':
                            # Only 1 word  is parsed after 'DELE ATOM', which
                            # is ok for "top_all27_prot_na.inp"
                            self.add(self.topology_info[cur_res][
                                     'DELE'], 'ATOM', words[2])
                        elif prop == 'ANGL':
                            # Need further parsing
                            self.add(self.topology_info[cur_res][
                                     'DELE'], 'ANGL', words[2:])
                        else:
                            # print 'WARNING!\n'+line+' NOT PARSED!\n' #
                            # Current "top_all27_prot_na.inp" file only has
                            # DELE ATOM and DELE ANGL, and therefore other
                            # situations are not parsed but warned
                            pass
        return error

    def patch_charmm_residue_atoms(self, residue, patch):
        '''
        patch 'patch' to 'residue' based on the charmm topology file
        '''
        # OPEN: ONLY ATOMS ARE PATCHED (NO BOND, ETC)
        error = []
        #
        # Make sure residue and patch are the right type
        #
        if residue not in self.topology_info:
            error.append('Residue ' + residue +
                         ' not found in the topology information list during patching!')
            return error
        elif patch not in self.topology_info:
            error.append(
                'Patch ' + residue + ' not found in the topology information list during patching!')
            return error
        #
        # Delete atoms from 'DELE' in patch list
        #
        tmp_dict = copy.deepcopy(self.topology_info[residue])
        # print 'ZHL ', patch, self.topology_info[patch]
        if 'DELE' in self.topology_info[patch]:  # ['DELE']:
            for atom_delete in self.topology_info[patch]['DELE']['ATOM']:
                for atom_residue in tmp_dict['ATOM']:
                    if atom_delete == atom_residue[0]:
                        tmp_dict['ATOM'].remove(atom_residue)
        #
        # Add atoms
        #
        num_added = 0
        for atom_patch in self.topology_info[patch]['ATOM']:
            for atom_residue in tmp_dict['ATOM']:
                if atom_residue[0] == atom_patch[0]:
                    tmp_dict['ATOM'].remove(atom_residue)
            if patch in ['NTER', 'GLYP', 'PROP']:
                tmp_dict['ATOM'].insert(num_added, atom_patch)
            if patch == 'CTER':
                tmp_dict['ATOM'].append(atom_patch)
            num_added += 1
        #
        # Add the patched residue to the topology dictionary as a new key of eg 'ALA_NTER'
        #
        self.topology_info[residue + '_' + patch] = tmp_dict
        self.charmm_residue_atoms[
            residue + '_' + patch] = numpy.array(tmp_dict['ATOM'])[:, 0].tolist()

    def setup_cys_patch_atoms_dirty(self):
        '''
        a dirty way to set up cys patch atoms
        simply remove HG1 atom from the atom list
        '''
        atoms = self.topology_info['CYS']['ATOM']
        my_atoms = []
        for atom in atoms:
            if atom[0] == 'HG1':
                pass
                # atoms.remove(atom)
            else:
                my_atoms.append(atom[0])

        return my_atoms

    def setup_charmm_residue_atoms(self):
        '''
        build the atom list of all the residues in the charmm topology file
        '''
        self.charmm_residue_atoms = {}
        for (key, value) in zip(self.topology_info.keys(), self.topology_info.values()):
            if type(value) is dict:
                if 'ATOM' in value:
                    atoms = numpy.array(value['ATOM'])[:, 0].tolist()
                    self.charmm_residue_atoms[key] = atoms
        # OPEN Setup CYS patch atoms in a dirty way

        self.charmm_residue_atoms['DISU'] = self.setup_cys_patch_atoms_dirty()

    def compare_list_ignore_order(self, l1, l2):
        '''
        compare two list ignoring the order
        '''
        if len(l1) != len(l2):  # make sure the lenght of two lists are the same
            return False
        for item in l1:
            # Dont allow duplicate elements during compare
            if l2.count(item) != 1:
                return False
        return True

    def check_charmm_atomic_order_reorganize(self):
        '''
        re-organize the atomic list according to the charmm topology contract
        do nothing if the atomic order alreay match that in the charmm topology file
        meanwhile make sure there are no missing or extra atoms
        H-atoms are required
        Patch N-ter for the first residue and C-ter for the last residue in each segment
        '''
        error = []
        bin_path = sasconfig.__bin_path__
        self.read_charmm_topology(topology_file_path=bin_path + '/toppar/')
        self.setup_charmm_residue_atoms()
        self.initialize_children()
        children_segname = self.init_child('segnames')
        for child_segname in children_segname:
            child_segname.initialize_children()
            children = child_segname.init_child('resids')
            resid_nter = child_segname.resid()[0]
            resid_cter = child_segname.resid()[-1]
            for child in children:
                child_resname = child.resname()[0]
                child_resid = child.resid()[0]
                child_names = child.name()
                child_indices = child.index()
                #
                # print child_resname,self.charmm_residue_atoms[child_resname], child_names
                # print
                # self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname],
                # child_names)
                if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                    if child_resname == 'CYS':
                        # OPEN: try DISU if CYS doesnt work
                        child_resname = 'DISU'
                    if child_resname == 'HIS':
                        # OPEN: try HSE if HIS doesnt work
                        child_resname = 'HSE'
                        if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                            # OPEN: try HSD if HIS doesnt work
                            child_resname = 'HSD'
                            if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                                # OPEN: try HSP if HIS doesnt work
                                child_resname = 'HSP'
                    if child_resid == resid_nter:
                            # print child_resname
                        if child_resname == 'GLY':
                            patch = 'GLYP'
                        elif child_resname == 'PRO':
                            patch = 'PROP'
                        elif child_resname == 'PROP':
                            patch = 'PROP'
                        elif child_resname == 'ADE' or child_resname == 'CTY' or child_resname == 'THY' or child_resname == 'GUA' or child_resname == 'URA':
                            pass

                        else:
                            patch = 'NTER'
                        self.patch_charmm_residue_atoms(child_resname, patch)
                        child_resname = child_resname + '_' + patch
                    if child_resid == resid_cter:
                        self.patch_charmm_residue_atoms(child_resname, 'CTER')
                        child_resname = child_resname + '_CTER'
                    if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                        error.append("For residue: " + child_resname + "\nthe atom name(s) do not match those in the charmm topology file!\n: charmm topology = " + str(
                            self.charmm_residue_atoms[child_resname]) + "\n your file = " + str(child_names))
                        return error
                #
                # if self.charmm_residue_atoms[child_resname] == child_names:
                #	continue
                # import pprint
                # pprint.pprint(child_resname)
                # pprint.pprint(self.charmm_residue_atoms[child_resname],width=100)
                new_indices = []
                for name in self.charmm_residue_atoms[child_resname]:
                    new_indices.append(child_names.index(name))
                if len(child_indices) != len(new_indices):
                    error.append('Number of atoms doesnt match that in charmm topology file for \nresname: ' +
                                 child_resname + '\nand resid: ' + child_resid + '!')
                    return error
                for i in range(len(child_indices)):
                    # if child_indices[i]!=new_indices[i]:
                    self.atom()[child_indices[i] -
                                1] = child.atom()[new_indices[i]]
                    # self.index()[indices[i]-1] = child.index[new_indices[i]]
                    self.name()[child_indices[i] -
                                1] = child.name()[new_indices[i]]
                    self.loc()[child_indices[i] -
                               1] = child.loc()[new_indices[i]]
                    self.resname()[child_indices[i] -
                                   1] = child.resname()[new_indices[i]]
                    self.chain()[child_indices[i] -
                                 1] = child.chain()[new_indices[i]]
                    self.resid()[child_indices[i] -
                                 1] = child.resid()[new_indices[i]]
                    self.rescode()[child_indices[i] -
                                   1] = child.rescode()[new_indices[i]]
                    self.occupancy()[child_indices[i] -
                                     1] = child.occupancy()[new_indices[i]]
                    self.beta()[child_indices[i] -
                                1] = child.beta()[new_indices[i]]
                    self.segname()[child_indices[i] -
                                   1] = child.segname()[new_indices[i]]
                    self.element()[child_indices[i] -
                                   1] = child.element()[new_indices[i]]
                    self.charge()[child_indices[i] -
                                  1] = child.charge()[new_indices[i]]
                    # only frame-0 was handled in sassubset.init_child
                    self.coor()[0][child_indices[i] -
                                   1] = child.coor()[0][new_indices[i]]
        return error

##### BEGIN ADDITIONS FROM ZAZMOL
##### BEGIN ADDITIONS FROM ZAZMOL
##### BEGIN ADDITIONS FROM ZAZMOL

    def renumber(self, **kwargs):
        '''
        Method to renumber index and resid fields

        Parameters
        ----------
        kwargs
        index :
        resid :

        Returns
        -------
        None
        updated system object
        index
        resid

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.index()[0]
        1
        >>> molecule.resid()[0]
        1

        default renumber() with no arguments renumbers
        both index and resid starting at 1

        >>> molecule.renumber()
        >>> molecule.index()[0]
        1
        >>> molecule.resid()[0]
        1

        passing a kwarg with a value will renumber appropriate
        attribute with value passed

        >>> molecule.renumber(index=3)
        >>> molecule.index()[0]
        3
        >>> molecule.resid()[0]
        1

        >>> molecule.renumber(resid=8)
        >>> molecule.index()[0]
        3
        >>> molecule.resid()[0]
        8

        >>> molecule.renumber(index=223, resid=18)
        >>> molecule.index()[0]
        223
        >>> molecule.resid()[0]
        18

        '''

        renumber_index_flag = False
        renumber_index_start = 1

        renumber_resid_flag = False
        renumber_resid_start = 1

        try:
            if kwargs['index']:
                renumber_index_flag = True
                renumber_index_start = kwargs['index']
        except:
            pass

        try:
            if kwargs['resid']:
                renumber_resid_flag = True
                renumber_resid_start = kwargs['resid']
        except:
            pass

        if renumber_index_flag:
            start = renumber_index_start
            index = numpy.array(
                [x for x in xrange(start, start + self._natoms)])
            self._index = index

        if renumber_resid_flag:
            resid = self._resid
            resid_array = []
            count = renumber_resid_start
            for i in xrange(len(resid)):
                this_resid = resid[i]
                if(i == 0):
                    last_resid = this_resid
                else:
                    if(this_resid != last_resid):
                        count += 1
                resid_array.append(count)
                last_resid = this_resid

            self._resid = numpy.array(resid_array, numpy.int)

        return

    def make_constraint_pdb(self, filename, basis_type, **kwargs):
        '''
     
        Method to rename attribute fields and assign values
        of 1.00 for the given basis type.

        "beta" is the default attribute field to reassign
        
        Default usage sets all attribute fields to zero
        prior to re-assigning basis type atoms to 1.00

        Parameters
        ----------
        filename    string : name of output PDB file to write
             
        basis types:

            'heavy : all atoms except hydrogen
          
            'protein' : all atoms in moltype protein
             
            'nucleic' : all atoms in moltype nucleic
                                        
            'backbone'
                -> proteins: N, CA, C, O
            
            'solute'    : all protein and nucleic


        kwargs
            optional arguments
                
                field='beta' : write 1.00 to beta field

                field='occupancy' : write 1.00 to occupancy field

                reset=True : set all values to 0.00 before setting
                                 individual vaules to 1.00

                reset=False: only change desired values to 1.00
                                 ignoring all other values
                    
        Returns
        -------
        None
            updated system object
            writes pdb to disk

        Examples
        --------
        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.make_constraint_pdb('constrainted_backbone_hiv1_gag.pdb', 'backbone')

        >>> molecule.make_constraint_pdb('constrainted_backbone_occupancy_hiv1_gag.pdb', 'backbone', field='occupancy')

        constrtain all heavy atoms

        >>> molecule.make_constraint_pdb('constrainted_heavy_hiv1_gag.pdb', 'heavy')

        constrain all protein atoms
        
        >>> molecule.make_constraint_pdb('constrainted_protein_hiv1_gag.pdb', 'protein')

        constrain all non-solvent atoms
        
        >>> molecule.make_constraint_pdb('constrainted_solute_hiv1_gag.pdb', 'solute')
        
        Note

        Only moltype `protein` and `nucleic` are supported.

        Output PDB file will contain coordinates from frame 0
             
        '''

        try:
            field = kwargs['field']
        except:
            field = 'beta'

        try:
            reset = kwargs['reset']
        except:
            reset = True
    
        if reset:
            new_value = ['0.00' for x in xrange(self._natoms)]
            if field == 'beta':
                self._beta = new_value
            elif field == 'occupancy':
                self._occupancy = new_value

        if basis_type == 'backbone':
            basis_string = "name[i] == 'N' or name[i] == 'CA' or name[i] == 'C' or name[i] == 'O'"
                           
        elif basis_type == 'heavy':
            basis_string = "name[i][0] != 'H'"

        elif basis_type == 'protein':
            basis_string = "moltype[i] == 'protein'"
        
        elif basis_type == 'nucleic':
            basis_string = "moltype[i] == 'nucleic'"

        elif basis_type == 'solute':
            basis_string = "moltype[i] == 'protein' or moltype[i] == 'nucleic'"

        error, mask = self.get_subset_mask(basis_string)

        if len(error) > 0:
            print('error = ', error)
            return

        if field == 'beta':
            descriptor = self._beta
        elif field == 'occupancy':
            descriptor = self._occupancy

        error = self.set_descriptor_using_mask(mask, descriptor, '1.00')
        
        if len(error) > 0:
            print('error = ', error)
            return
               
        if field == 'beta':
            self._beta = descriptor
        elif field == 'occupancy':
            self._occupancy = descriptor

        self.write_pdb(filename, 0, 'w')

        return


    def make_backbone_pdb_from_fasta(self, filename, moltype, **kwargs):
        """
        Method to write a PDB file of one atom per residue
        based on input FASTA sequence

        http://en.wikipedia.org/wiki/FASTA_format

        Parameters
        ----------
        filename    string : name of output PDB file to write
        
        moltype string
                -> 'protein' ; fasta sequence of a protein
                -> 'nucleic' ; fasta sequence of a nucleic acid
             
        kwargs
            optional future arguments
                                                                                     
        Returns
        -------
        self._fasta
            sets the fasta attribute in a system object

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule("hiv1_gag.pdb")
        >>> molecule.create_fasta()
        >>> molecule.fasta()[:5]
        ['G', 'A', 'R', 'A', 'S']

        >>> molecule.make_backbone_pdb_from_fasta('sequence.pdb', 'protein')

        Note
        ____

        Only protein and nucleic acid one-letter codes are supported
        Protein PDB files are written with one "CA" atom per residue
        Nucleic acid PDF files are written with on "O5'" atom per residue
       
        All coordinate values are set to 0.000
         
        """

        if moltype == "protein":
            residue_dictionary = self.one_to_three_letter_protein_residue_dictionary
            sequence_name = "CA"

        elif moltype == "nucleic":
            residue_dictionary = self.one_to_three_letter_nucleic_residue_dictionary
            sequence_name = "O5'"

        sequence = []
        first_flag = True

        for residue in self._fasta:

            if moltype == "protein":
                if residue == 'G' and first_flag:
                    #sequence.append('GLYP')
                    sequence.append('GLY')
                elif residue == 'P' and first_flag:
                    #sequence.append('PROP')
                    sequence.append('PRO')
                else:
                    sequence.append(residue_dictionary[residue])

            elif moltype == "nucleic":
                sequence.append(residue_dictionary[residue])
        
            first_flag = False

        #import sasmol.system as system
        #molecule = system.Molecule_Maker(len(sequence), name=sequence_name)

        #import sasmol.sasmol as sasmol
        molecule = sasmol.Molecule_Maker(len(sequence), name=sequence_name)

        molecule._resname = sequence

        molecule.write_pdb(filename, 0, 'w')

        return

    def create_fasta(self, **kwargs):
        """
        Method to make a fasta file compatible lists of residue names

        http://en.wikipedia.org/wiki/FASTA_format

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        self._fasta
            sets the fasta attribute in a system object 

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule("hiv1_gag.pdb")
        >>> molecule.create_fasta()
        >>> print(molecule.fasta()[:5]) # doctest: +NORMALIZE_WHITESPACE
        ['G', 'A', 'R', 'A', 'S'] 

        >>> molecule.create_fasta(fasta_format=True)
        >>> print(molecule.fasta()[:5])
        >
        GAR
        >>> molecule.create_fasta(fasta_format=True,width='60')
        >>> print(molecule.fasta()[:-1])
        >
        GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQIL
        GQLQPSLQTGSEELRSLYNTIAVLYCVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAAD
        TGNNSQVSQNYPIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATP
        QDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTS
        TLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFY
        KTLRAEQASQEVKNAATETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKAR
        VIAEAMSQVTNSATIMMQKGNFRNQRKTVKCFNCGKEGHIAKNCRAPRKKGCWKCGKEGH
        QMKDCTERQAN

        >>> molecule.create_fasta(fasta_format=True,width='60',name='aar')
        >>> print(molecule.fasta()[:90])
        >aar
        GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQIL
        GQLQPSLQTGSEELRSLYNTIAVL

        
        Note
        ----
        Other use types below

        molecule.create_fasta(fasta_format=True,exclude_hetatm=True,by_chain=True)
        print molecule.fasta()

        print '>>> testing by_chain with HETATM (t.py): '

        molecule.create_fasta(fasta_format=True,by_chain=True)
        print molecule.fasta()

        print '>>> testing by_segname (t.py): '

        molecule.create_fasta(fasta_format=True,exclude_hetatm=True,by_segname=True)
        print molecule.fasta()

        Note that this creates a simple string that is associated with the molecule (self)
        and it will return without assigning a string if a non-standard three letter code
        is supplied.

        """

        fasta_format = False
        by_chain = False
        by_segname = False
        single_chain = False
        single_segname = False
        exclude_hetatm = False

        if 'fasta_format' in kwargs:
            fasta_format = True

        if 'name' in kwargs:
            name = kwargs['name']
            header = '>' + name
        else:
            header = '>'

        if 'width' in kwargs:
            width = kwargs['width']
        else:
            width = '80'

        if 'exclude_hetatm' in kwargs:
            exclude_hetatm = True

        if 'by_chain' in kwargs:
            by_chain = True

        elif 'by_segname' in kwargs:
            by_segname = True

        one_resname = []

        resname = self.resname()
        atom = self.atom()
        chain = self.chain()
        segname = self.segname()

        residue_dictionary = {
            'ALA': 'A',
            'ARG': 'R',
            'ASP': 'D',
            'ASN': 'N',
            'CYS': 'C',
            'GLU': 'E',
            'GLN': 'Q',
            'GLY': 'G',
            'HSD': 'H',
            'HIS': 'H',
            'HSE': 'H',
            'HSP': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V',
            'GUA': 'G',
            'CYT': 'C',
            'ADE': 'A',
            'THY': 'T',
            'URA': 'U'
        }
 
        for i in xrange(len(resname)):
            this_resname = resname[i]
            if this_resname in residue_dictionary:
                one_resname.append(residue_dictionary[this_resname])
            elif (atom[i] == 'HETATM'):
                # print 'skipping non-standard resname in HETATM record:
                # ',this_resname
                one_resname.append("X")
            else:
                print('non standard resname: ', this_resname)
                print('unable to make Fasta conversion')
                return

        self.setOne_letter_resname(one_resname)

        resid = self.resid()

        last_resid = None
        last_chain = None
        last_segname = None
        fasta = []
        local_fasta = []
        number_of_chains = 0
        number_of_segnames = 0
        chain_name = []
        segname_name = []
        first = True
        for i in xrange(len(resid)):
            this_resid = resid[i]

            if by_chain:
                this_chain = chain[i]
                if this_chain != last_chain:
                    number_of_chains += 1
                    chain_name.append(this_chain)
                    if first:
                        first = False
                    else:
                        fasta.append(local_fasta)

                    last_chain = this_chain
                    local_fasta = []

            if by_segname:
                this_segname = segname[i]
                if this_segname != last_segname:
                    number_of_segnames += 1
                    segname_name.append(this_segname)
                    if first:
                        first = False
                    else:
                        fasta.append(local_fasta)

                    last_segname = this_segname
                    local_fasta = []

            this_resname = self.one_letter_resname()[i]
            if this_resid != last_resid:
                local_fasta.append(this_resname)
                last_resid = this_resid

        if by_chain or by_segname:
            fasta.append(local_fasta)

        elif not by_chain and not by_segname:
            fasta = local_fasta
            number_of_chains = 1
            chain_name = chain[0]

        final_fasta = ''

        if by_segname:
            number_of_chains = number_of_segnames

        for i in xrange(number_of_chains):
            saveme = False
            if fasta_format:
                from textwrap import TextWrapper
                wrapper = TextWrapper(width=int(width))
                if by_chain or by_segname:
                    if exclude_hetatm:
                        while "X" in fasta[i]:
                            fasta[i].remove("X")

                    joined_fasta = ''.join(fasta[i])
                else:
                    if exclude_hetatm:
                        while "X" in fasta:
                            fasta.remove("X")
                    joined_fasta = ''.join(fasta)

                formatted_fasta = "\n".join(wrapper.wrap(joined_fasta))

                if(len(formatted_fasta.strip()) > 0):
                    saveme = True

                if(len(header) > 1):
                    if by_chain:
                        formatted_fasta = header + ' chain:' + \
                            chain_name[i] + '\n' + formatted_fasta
                    elif by_segname:
                        formatted_fasta = header + ' segname:' + \
                            segname_name[i] + '\n' + formatted_fasta
                    else:
                        formatted_fasta = header + '\n' + formatted_fasta

                else:
                    if by_chain:
                        formatted_fasta = header + 'chain:' + \
                            chain_name[i] + '\n' + formatted_fasta
                    elif by_segname:
                        formatted_fasta = header + 'segname:' + \
                            segname_name[i] + '\n' + formatted_fasta
                    else:
                        formatted_fasta = header + '\n' + formatted_fasta

                if saveme:
                    final_fasta += formatted_fasta + '\n'
                else:
                    final_fasta += '\n'

        if fasta_format:
            self.setFasta(final_fasta)
        else:
            self.setFasta(fasta)

        return

##### END ADDITIONS FROM ZAZMOL
##### END ADDITIONS FROM ZAZMOL
##### END ADDITIONS FROM ZAZMOL


def set_pdb_values(m, natoms, **kwargs):

    atom = []
    index = []
    name = []
    loc = []
    resname = []
    chain = []
    resid = []
    rescode = []
    occupancy = []
    beta = []
    segname = []
    element = []
    charge = []
    moltype = []

    if 'name' in kwargs:
        default_name = kwargs['name']
    else:
        default_name = "C"

    if 'resname' in kwargs:
        default_resname = kwargs['resname']
    else:
        default_resname = "DUM"

    if 'segname' in kwargs:
        default_segname = kwargs['segname']
    else:
        default_segname = "DUM"

    if 'element' in kwargs:
        default_element = kwargs['element']
    else:
        default_element = "C"

    if 'moltype' in kwargs:
        default_moltype = kwargs['moltype']
    else:
        default_moltype = "other"

    if 'coor' in kwargs:
        m.setCoor = kwargs['coor']

    for i in xrange(natoms):
        atom.append("ATOM  ")
        index.append(i + 1)
        name.append(default_name)
        loc.append(" ")
        resname.append(default_resname)
        chain.append(" ")
        resid.append(i + 1)
        rescode.append(" ")
        occupancy.append("  0.00")
        beta.append("  0.00")
        segname.append(default_segname)
        element.append(default_element)
        charge.append("  ")
        moltype.append(default_moltype)

    m.setAtom(atom)
    m.setIndex(index)
    m.setName(name)
    m.setLoc(loc)
    m.setResname(resname)
    m.setChain(chain)
    m.setResid(resid)
    m.setRescode(rescode)
    m.setOccupancy(occupancy)
    m.setBeta(beta)
    m.setSegname(segname)
    m.setElement(element)
    m.setCharge(charge)
    m.setMoltype(moltype)
    m.setNatoms(natoms)

    return
