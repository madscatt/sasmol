#include <sasmol.h>
#include <sasio.h>
#include <sascalc.h>
#include <sasop.h>
#include <sassubset.h>
#include <dcdio.h>
#include <util.h>
#include <vector>
#include <boost/dynamic_bitset.hpp>

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

	float rg = mol.calc_rg(frame) ;
	std::cout << "radius of gyration = " << rg << std::endl ;
	std::cout << "sasmol.py : radius of gyration = 64.043168998442439 " << std::endl ;

	mol.center(frame) ;
	
	mol.move_to(frame,com) ;

	mol.calc_pmi(frame) ;

	std::cout << "mol.uk = " << mol.uk << std::endl ;
	std::cout << "mol.ak = " << mol.ak << std::endl ;

	float cutoff = 0.8 ;
	int check = 1 ;  // 1 == overlap and 0 == no overlap

	std::cout << ">>> checking for overlap (xyz): cutoff = "<< cutoff << std::endl ;

	clock_t tStart = clock();

	check = mol.self_overlap(cutoff,frame) ;

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	if(check == 0) 
	{
		std::cout << ">>> NO OVERLAP FOUND <<<" << std::endl;
	}
	else
	{
		std::cout << ">>> OVERLAP FOUND <<<" << std::endl;
	}


	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;

	mol.write_dcd(dcd_output_filename) ;	

	mol.write_pdb(output_filename,frame);

	const std::string second_dcd_output_filename="second_new_min3.dcd" ;
	
	FILE *outfile = sasio::open_write_dcd_file(second_dcd_output_filename, mol.natoms, mol.number_of_frames);

	std::vector<float> value = {0.0,0.0,0.0} ;
	float vcount = 0.0 ;
	std::cout << value[0] << " " << value[1] << " " << value[2] << std::endl ;
	
	tStart = clock();
	int step = 0 ;

	const std::string translated_pdb_filename = "t_min3.pdb" ;

	for(int j = 0 ; j < 1000 ; ++j)
	{
		mol.translate(frame,value) ;	

		mol.write_dcd_step(outfile,step,j);

		value[0] += 0.0001 ;
		value[1] += 0.0001 ;
		value[2] += 0.0001 ;
	
		vcount += value[0] ;
		if(j == 999) mol.write_pdb(translated_pdb_filename,frame) ;
	}

	int close_result = close_dcd_write(outfile) ;
	
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	std::cout << vcount << std::endl;

	std::cout << value[0] << " " << value[1] << " " << value[2] << std::endl ;

	std::cout << "TESTING DCD READ" << std::endl << std::endl ;

	std::string dcd_input_file_name = "new_min3.dcd" ;
	int input_natoms, input_nset, input_reverseEndian ;

	FILE *dcd_file_pointer ;

	dcd_file_pointer = sasio::open_dcd_read(dcd_input_file_name, input_natoms, input_nset, input_reverseEndian) ;

	std::cout << "BACK FROM OPEN DCD READ" << std::endl ;

	std::cout << "input_natoms = " << input_natoms << std::endl ;
	std::cout << "input_nset = " << input_nset << std::endl ;
	std::cout << "input_reverseEndian = " << input_reverseEndian << std::endl ;

	std::cout << "TESTING DCD READ STEP" << std::endl << std::endl ;

	sasmol::SasMol mol2 ; 
	mol2.read_pdb(pdb_filename);
	int input_frame = 0 ;

	mol2.read_dcd_step( dcd_file_pointer, input_frame, input_natoms, input_reverseEndian) ;	
	
	int close_read_result = close_dcd_read(dcd_file_pointer) ;

	std::string second_pdb_filename = "mol2_new_min3.pdb" ;

	mol2.write_pdb(second_pdb_filename,frame);

	std::cout << "BACK FROM OPEN DCD READ" << std::endl ;

	std::cout << "TESTING DCD READ (FULL)" << std::endl << std::endl ;

	std::string second_dcd_filename = "second_new_min3.dcd" ;
	
	mol.read_dcd(second_dcd_filename);
	
	std::cout << "BACK FROM DCD READ (FULL)" << std::endl ;

	std::cout << "number_of_frames = " << mol.number_of_frames << std::endl ;


	std::cout << "checking RMSD calculator " << std::endl ;


	sasmol::SasMol mol3 ; sasmol::SasMol mol5 ;
	std::cout << "reading translated_pdb_filename" << std::endl ;
	mol3.read_pdb(translated_pdb_filename) ;
	std::cout << "reading pdb_filename" << std::endl ;
	mol5.read_pdb(pdb_filename) ;
	frame = 0 ;
	std::cout << "calculating rmsd" << std::endl ;
	float rmsd = mol3.calc_rmsd(mol5,frame) ;

	std::cout << "rmsd = " << rmsd << std::endl ;
	std::cout << "sasmol.py : rmsd = 121.09128746800168" << std::endl ;

	std::cout << "checking minmax calculator " << std::endl ;

	std::vector<std::vector<float >> minmax = mol5.calc_minmax(frame) ;	

	std::cout << "min values" << std::endl ; util::print(minmax[0]) ;
	std::cout << "max values" << std::endl ; util::print(minmax[1]) ;
	
	std::cout << "sasmol.py : minmax = [array([-31.29899979, -93.23899841, -85.81900024]), array([ 19.64699936,  30.37800026,  99.52999878])]" << std::endl ;

	mol.number_of_frames = 1 ;

	sasmol::SasMol mol4 ;

	mol4.read_pdb(pdb_filename);
	
	frame = 0 ;	
	std::cout << "mol4.x(0,0) = " << mol4.x(0,0) << std::endl ;
	std::cout << "mol4.x(6729,0) = " << mol4.x(6729,0) << std::endl ;

	com = mol4.calc_com(frame) ; std::cout << "frame 0:" << std::endl ; util::print(com) ;
	std::cout << "nf = " << mol4.number_of_frames << std::endl ;
	std::cout << "reading pdb" << std::endl;	
	mol4.read_pdb(translated_pdb_filename);
	std::cout << "nf = " << mol4.number_of_frames << std::endl ;

	std::cout << "mol4.x(0,0) = " << mol4.x(0,0) << std::endl ;
	std::cout << "mol4.x(6729,0) = " << mol4.x(6729,0) << std::endl ;

	std::cout << "mol4.x(0,1) = " << mol4.x(0,1) << std::endl ;
	std::cout << "mol4.x(6729,1) = " << mol4.x(6729,1) << std::endl ;

	com = mol4.calc_com(frame) ; std::cout << "frame " << frame << " again:" << std::endl ;  util::print(com) ;
	frame++ ;
	com = mol4.calc_com(frame) ; std::cout << "frame " << frame << " :" << std::endl ;  util::print(com) ;
	
	std::cout << std::endl << "testing find_u" << std::endl << std::endl ;

	frame = 0 ;

	sasmol::SasMol mol6 ;
	mol6.read_pdb(pdb_filename) ;

	Eigen::Matrix3f u = mol6.find_u(mol3,frame) ;

	std::cout << ">>> U = \n" << u << std::endl << std::endl ;
	std::cout << "sasmath.py : u = \n-0.945133  0.310686  0.100992\n-0.300624 -0.706126 -0.641102\n-0.127869 -0.636287  0.760782 \n" << std::endl ;

	std::cout << std::endl << "testing align" << std::endl << std::endl ;

	sasmol::SasMol mol10 ; sasmol::SasMol mol11 ; mol10.read_pdb(pdb_filename) ; mol11.read_pdb(translated_pdb_filename) ;

	mol10.align(mol11,frame) ;	// align min3.pdb ONTO t_min3.pdb

	const std::string aligned_10_filename = "aligned_10_min3.pdb" ;
	mol10.write_pdb(aligned_10_filename,frame);

	const std::string aligned_11_filename = "aligned_11_min3.pdb" ;
	mol11.write_pdb(aligned_11_filename,frame);
	
	std::cout << std::endl << "testing signed angle" << std::endl << std::endl ;
	
	sasmol::SasMol mol7 ;
	mol7.read_pdb(pdb_filename) ;

	Eigen::Matrix<float,3,1> coor1 ; coor1 << mol7.x(0,0),mol7.y(0,0),mol7.z(0,0) ;
	Eigen::Matrix<float,3,1> coor2 ; coor2 << mol7.x(1,0),mol7.y(1,0),mol7.z(1,0) ;
	Eigen::Matrix<float,3,1> coor3 ; coor3 << mol7.x(2,0),mol7.y(2,0),mol7.z(2,0) ;
	Eigen::Matrix<float,3,1> coor4 ; coor4 << mol7.x(3,0),mol7.y(3,0),mol7.z(3,0) ;

	float sa = sasmath::signed_angle(coor1,coor2,coor3) ;

	std::cout << ">>> signed_angle = \n" << sa << std::endl << std::endl ;
	std::cout << "sasmath.py : sa = \n-0.363156\n" << std::endl ;

	std::cout << std::endl << "testing signed angle" << std::endl << std::endl ;
	
	float da = sasmath::dihedral_angle(coor1,coor2,coor3,coor4) ;

	std::cout << ">>> dihedral angle = \n" << da << std::endl << std::endl ;
	std::cout << "sasmath.py : da = \n-35.8585807936\n" << std::endl ;

	std::cout << std::endl << "testing calc_angle" << std::endl << std::endl ;
	
	float angle = sasmath::calc_angle(coor1,coor2,coor3) ;

	std::cout << ">>> calc_angle = \n" << angle << std::endl << std::endl ;

	std::cout << "sasmath.py : angle = \n0.667242537186\n" << std::endl ;


	sasmol::SasMol mol8 ;
	mol8.read_pdb(translated_pdb_filename);

	std::cout << std::endl << "testing rotate" << std::endl << std::endl ;

	std::string axis = "x" ;
	float theta = acos(-1.0)/4.0 ;
	mol7.rotate(frame,axis,theta) ;
	
	const std::string rotated_filename = "rotated_min3.pdb" ;
	mol7.write_pdb(rotated_filename,frame);

	std::cout << "atom 1(xyz) = " << mol7.x(0,0) << " " << mol7.y(0.0) << " " << mol7.z(0.0) << std::endl ;
	std::cout << "atom 2(xyz) = " << mol7.x(1,0) << " " << mol7.y(1.0) << " " << mol7.z(1.0) << std::endl ;
	std::cout << "atom 3(xyz) = " << mol7.x(2,0) << " " << mol7.y(2.0) << " " << mol7.z(2.0) << std::endl ;

	std::cout << "SHOULD BE : \n" ;

	std::cout << "atom 1(xyz) = -21.525 -110.364 32.7058\n" ;
	std::cout << "atom 2(xyz) = -22.003 -111.093 32.1649\n" ;
	std::cout << "atom 3(xyz) = -21.905 -110.458 33.695\n\n" ;
	std::cout << " ... and w/o com correction ..." << std::endl ;
	std::cout << "atom 1 (xyz) -21.525-109.121  13.574\n" ;
	std::cout << "atom 2 (xyz) -22.003-109.850  13.033\n" ;
	std::cout << "atom 3 (xyz) -21.905-109.215  14.564\n\n" ;

	util::pp("now work on sassubset") ;




	const std::string test_filename = "all_test.pdb" ;
	sasmol::SasMol mol22 ;

	mol22.read_pdb(test_filename) ;

	mol22.write_pdb("did_it_work.pdb",frame) ;


	std::cout << std::endl << std::endl ;
	util::print_run_details() ;
	std::cout << std::endl << " >>> DONE <<< " << std::endl << std::endl ;
	return 0 ;


}
