#include <sasmol.h>
#include <sasio.h>
#include <sascalc.h>
#include <sasop.h>
#include <SasPteros/saspteros.h>
#include <dcdio.h>
#include <util.h>
#include <algorithm>
#include <iostream>

#include <vector>

/********* methods        ******************/

/********* main           ******************/

int main(){

	std::cout << "\n\n\n" ;
	
	sasmol::SasMol mol ;
	
	util::compile_time(__FILE__, __func__, __DATE__, __TIME__ );

//	const std::string filename="ubq.xyz" ;
//	mol.read_xyz(filename) ;
      
	const std::string pdb_filename = "min3.pdb" ;
	const std::string output_filename = "new_min3.pdb" ;
	const std::string dcd_output_filename="new_min3.dcd" ;

	std::cout << "testing sasmol " << " " << pdb_filename << std::endl ; 

	mol.read_pdb(pdb_filename);

	int frame = 0 ;

	std::cout << "number of atoms = " << mol._natoms() << std::endl ;
	std::cout << "total mass = " << mol._total_mass() << std::endl ;
	
	std::vector<float> com = mol.calc_com(frame) ;
	util::print(com) ;	

    /// test for get_subset_mask
    std::cout<<std::endl<<"Testing get_subset_mask..."<<std::endl;
    boost::dynamic_bitset<> mask = mol.get_subset_mask(std::string("resid 3"));
    for (boost::dynamic_bitset<>::size_type i = 0; i<std::min(int(mask.size()),100); ++i) std::cout<<mask[i]<<" "; std::cout<<std::endl;

    /// test for get_coor_using_mask
    std::cout<<std::endl<<"Testing get_coor_using_mask..."<<std::endl;
    Eigen::Matrix<float,3,Eigen::Dynamic> coor = mol.get_coor_using_mask(0, mask);
    std::cout<<coor<<std::endl;

	sasmol::SasMol mol2;
	mol2.read_pdb("t_min3.pdb");

    /// test for set_coor_using_mask
    std::cout<<std::endl<<"Testing set_coor_using_mask..."<<std::endl;
    mol.set_coor_using_mask(mol2, 0, mask);
    coor = mol.get_coor_using_mask(0, mask);
    std::cout<<"After set_coor_using_mask: "<<std::endl;
    std::cout<<coor<<std::endl;

    /// test for copy_molecule_using_mask
    std::cout<<"Testing copy_molecule_using_mask..."<<std::endl;
    sasmol::SasMol mol3;
    mol.copy_molecule_using_mask(mol3, mask, 0);
    mol3.write_pdb("min3_copy_molecule_using_mask.pdb",frame);
    std::cout<<"results written to min3_copy_molecule_using_mask.pdb"<<std::endl;

    /// test for merge_two_molecules
    std::cout<<std::endl<<"Testing merge_two_molecules..."<<std::endl;
    sasmol::SasMol mol4;
    mol4.read_pdb("min3_copy_molecule_using_mask.pdb");
    sasmol::SasMol mol5;
    mol5.merge_two_molecules(mol, mol4);
    mol5.write_pdb("min3_merge_two_molecules.pdb",frame);
    std::cout<<"results written to min3_merge_two_molecules.pdb"<<std::endl;

    /// test for get_dihedral_subset_mask
    std::cout<<std::endl<<"Testing get_dihedral_subset_mask..."<<std::endl;
    std::vector<int> flexible_residues={2,5};
    std::vector<boost::dynamic_bitset<>> mask_dihedral = mol.get_dihedral_subset_mask(flexible_residues, 0);
    for (int i=0; i<mask_dihedral.size(); ++i) {for (int j=0; j<std::min(int(mol._natoms()),100); ++j) std::cout<<mask_dihedral[i][j]<<" "; std::cout<<std::endl;}

    /// test for initialize_children
    std::cout<<std::endl<<"Testing initialize_children..."<<std::endl;
    mol.initialize_children();
    for (int i=0; i<mol._mask_resnames().size(); ++i) {std::cout<<mol._unique_resname()[i].c_str()<<" "; for (int j=0; j<std::min(int(mol._natoms()),100); ++j) std::cout<<mol._mask_resnames()[i][j]<<" "; std::cout<<std::endl;}

    /// test for duplicate_molecule
    std::cout<<std::endl<<"Testing duplicate_molecule..."<<std::endl;
    sasmol::SasMol mol6;
    com = {0.,0.,0.};
    const int frame_duplicate = 10;
    std::vector<std::vector<float>> com_coor;
    for (int i=0; i<frame_duplicate; ++i) com_coor.push_back(com);
    bool flag_rotate = true;
    mol.duplicate_molecule(mol6, frame_duplicate, frame, com_coor, flag_rotate);
    int frame_buggy = 0;
    mol6.write_pdb("min3_duplicate_molecule.pdb",frame_buggy);
    std::cout<<"results written to min3_duplicate_molecule.pdb"<<std::endl;

    /// test for set_descriptor_using_mask
    std::cout<<std::endl<<"Testing set_descriptor_using_mask..."<<std::endl;
    mol.set_descriptor_using_mask(mask, mol._atom_segname(), std::string("T"));
    mol.write_pdb("min3_set_descriptor_using_mask.pdb",frame);
    std::cout<<"results written to min3_set_descriptor_using_mask.pdb"<<std::endl;

	std::cout << std::endl << std::endl ;
	util::print_run_details() ;
	std::cout << std::endl << " >>> DONE <<< " << std::endl << std::endl ;
    return 0;
}
