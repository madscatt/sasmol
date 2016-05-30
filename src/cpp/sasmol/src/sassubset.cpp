#include <dcdio.h>
#include <sasio.h>
#include <sasmol.h>
#include <sascalc.h>
#include <sassubset.h>
#include <SasPteros/saspteros.h>
#include <util.h>
#include <boost/dynamic_bitset.hpp>
#include <stdexcept> 
#include <stdlib.h>
#include <time.h>

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void 
sassubset::Mask::
get_mask_array(std::vector<boost::dynamic_bitset<>> & farray, int nflexible, int natoms, const std::vector<std::string> & name, const std::vector<int> & resid, const std::vector<int> & flexible_residues, const int nresidues, const int mtype)
{

	int i, j, q0, fr, value, count ;
	int value_array[natoms] ;

	for(i=0;i<natoms;i++){
		value_array[i]=0 ;
	}

        for(i=0;i<nflexible;i++) {
		for(j=0;j<natoms;j++){
			farray[i][j]=0;
		}
	}

	count=0 ;

	if(mtype == 0) {
		for(fr=0 ; fr<nflexible ; fr++){
			q0 = flexible_residues[fr] ;
			for(i=0 ; i<natoms ; i++){
				value=0 ;
				if (resid[i]<q0-1 || resid[i]>q0+1)
                		{
                    			count++;
                    			value_array[i]=0;
                    			continue;
                		}
				if(resid[i] == q0-1 && (strcmp(name[i].c_str(),"C")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"N")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"CA")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"C")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i].c_str(),"N")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else{
					value_array[i]=value ;
				}
			}
			for(i=0 ; i<natoms ; i++){
				farray[fr][i]=value_array[i] ;		
			} 
		} 
	}
	else if(mtype == 1){
		for(fr=0 ; fr<nflexible ; fr++){
			q0 = flexible_residues[fr] ;
			for(i=0 ; i<natoms ; i++){
				value=0 ;
				if (resid[i]<q0-1 || resid[i]>q0+1)
                		{
                    			count++;
                    			value_array[i]=0;
                    			continue;
                		}
				if(resid[i] == q0-1 && (strcmp(name[i].c_str(),"O3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"P")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"O5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"C5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"C4'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"C3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i].c_str(),"O3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i].c_str(),"P")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i].c_str(),"O5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else{
					value_array[i]=value ;
				}
			}
			for(i=0 ; i<natoms ; i++){
				farray[fr][i]=value_array[i] ;		
			} 
		} 
	}

	return ;

}
	
///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
sassubset::Mask &
sassubset::Mask::
initialize_children()
{
    const int number_of_names = _unique_name().size();
    const int number_of_resnames = _unique_resname().size();
    const int number_of_resids = _unique_resid().size();
    const int number_of_chains = _unique_chain().size();
    const int number_of_occupancies = _unique_occupancy().size();
    const int number_of_betas = _unique_beta().size();
    const int number_of_elements = _unique_element().size();
    const int number_of_segnames = _unique_segname().size();

    const int natoms = _natoms();
    int i;
    for (i=0; i<number_of_names; ++i) mask_names.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_resnames; ++i) mask_resnames.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_resids; ++i) mask_resids.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_chains; ++i) mask_chains.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_occupancies; ++i) mask_occupancies.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_betas; ++i) mask_betas.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_elements; ++i) mask_elements.push_back(boost::dynamic_bitset<>(natoms));
    for (i=0; i<number_of_segnames; ++i) mask_segnames.push_back(boost::dynamic_bitset<>(natoms));

    for (i=0; i<natoms; ++i)
    {
        mask_names[std::find(_unique_name().begin(),_unique_name().end(),_atom_name()[i])-_unique_name().begin()][i] = 1;
        mask_resnames[std::find(_unique_resname().begin(),_unique_resname().end(),_atom_resname()[i])-_unique_resname().begin()][i] = 1;
        mask_resids[std::find(_unique_resid().begin(),_unique_resid().end(),_atom_resid()[i])-_unique_resid().begin()][i] = 1;
        mask_chains[std::find(_unique_chain().begin(),_unique_chain().end(),_atom_chain()[i])-_unique_chain().begin()][i] = 1;
        mask_occupancies[std::find(_unique_occupancy().begin(),_unique_occupancy().end(),_atom_occupancy()[i])-_unique_occupancy().begin()][i] = 1;
        mask_betas[std::find(_unique_beta().begin(),_unique_beta().end(),_atom_beta()[i])-_unique_beta().begin()][i] = 1;
        mask_elements[std::find(_unique_element().begin(),_unique_element().end(),_atom_element()[i])-_unique_element().begin()][i] = 1;
        mask_segnames[std::find(_unique_segname().begin(),_unique_segname().end(),_atom_segname()[i])-_unique_segname().begin()][i] = 1;
    }

    return * this;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<boost::dynamic_bitset<>>
sassubset::Mask::
get_dihedral_subset_mask(const std::vector<int> & flexible_residues, const int mtype)
{
    const int natoms = _natoms();
	const int nflexible = flexible_residues.size();
    const int nresidues = _atom_resid()[natoms-1] - _atom_resid()[0] + 1;
    std::vector<boost::dynamic_bitset<>> masks;
    for (int i = 0; i < nflexible; ++i) masks.push_back(boost::dynamic_bitset<>(natoms));
    std::vector<std::string> name = util::strip_white_space(_atom_name());
	get_mask_array(masks, nflexible, natoms, name, _atom_resid(), flexible_residues, nresidues, mtype);
	return masks;

}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
std::vector<int>
sassubset::Mask::
get_indices_from_mask(const boost::dynamic_bitset<> & mask)
{
    std::vector<int> index;
    for (boost::dynamic_bitset<>::size_type i = 0; i<mask.size(); ++i)
    {
        if(mask[i]) index.push_back(i);
    }
    return index;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
boost::dynamic_bitset<>
sassubset::Mask::
get_mask_from_indices(const int natoms, const std::vector<int> & index)
{
    boost::dynamic_bitset<> mask(natoms);
    for (std::vector<int>::size_type i=0; i<index.size(); ++i)
    {
        mask[index[i]] = 1;
    }
    return mask;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
boost::dynamic_bitset<>
sassubset::Mask::
get_subset_mask(const std::string & basis_filter)
{
    int natoms = _natoms();

    SasPteros sel(dynamic_cast<sasmol::SasMol &>(*this), basis_filter);
    std::vector<int> index = sel.index;
    boost::dynamic_bitset<> mask = get_mask_from_indices(natoms, index);

    return mask;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sassubset::Mask::
merge_two_molecules(sasmol::SasMol & mol1, sasmol::SasMol & mol2)
{
    _atom_record().clear();
    _atom_index().clear();
    _atom_name().clear();
    _atom_altloc().clear();
    _atom_resname().clear();
    _atom_chain().clear();
    _atom_resid().clear();
    _atom_icode().clear();
    _atom_occupancy().clear();
    _atom_beta().clear();
    _atom_segname().clear();
    _atom_element().clear();
    _atom_selement().clear();
    _atom_charge().clear();
    /// @note to ZHL: don't know how to clear coor yet
    for (int i = 0; i < mol1._natoms(); ++i)
    {
        _atom_record().push_back(mol1._atom_record()[i]);
        _atom_index().push_back(mol1._atom_index()[i]);
        _atom_name().push_back(mol1._atom_name()[i]);
        _atom_altloc().push_back(mol1._atom_altloc()[i]);
        _atom_resname().push_back(mol1._atom_resname()[i]);
        _atom_chain().push_back(mol1._atom_chain()[i]);
        _atom_resid().push_back(mol1._atom_resid()[i]);
        _atom_icode().push_back(mol1._atom_icode()[i]);
        _atom_occupancy().push_back(mol1._atom_occupancy()[i]);
        _atom_beta().push_back(mol1._atom_beta()[i]);    
        _atom_segname().push_back(mol1._atom_segname()[i]);    
        _atom_element().push_back(mol1._atom_element()[i]);    
        _atom_selement().push_back(mol1._atom_selement()[i]);    
        _atom_charge().push_back(mol1._atom_charge()[i]);    
    }
    const int last_index_mol1 = mol1._atom_index()[mol1._natoms()-1];
    int this_index = last_index_mol1+1;
    for (int i = 0; i < mol2._natoms(); ++i)
    {
        _atom_record().push_back(mol2._atom_record()[i]);
        _atom_index().push_back(this_index);
        _atom_name().push_back(mol2._atom_name()[i]);
        _atom_altloc().push_back(mol2._atom_altloc()[i]);
        _atom_resname().push_back(mol2._atom_resname()[i]);
        _atom_chain().push_back(mol2._atom_chain()[i]);
        _atom_resid().push_back(mol2._atom_resid()[i]);
        _atom_icode().push_back(mol2._atom_icode()[i]);
        _atom_occupancy().push_back(mol2._atom_occupancy()[i]);
        _atom_beta().push_back(mol2._atom_beta()[i]);    
        _atom_segname().push_back(mol2._atom_segname()[i]);    
        _atom_element().push_back(mol2._atom_element()[i]);    
        _atom_selement().push_back(mol2._atom_selement()[i]);    
        _atom_charge().push_back(mol2._atom_charge()[i]);    
        ++ this_index;
    }
    _natoms() = _atom_index().size();
    const int frame = 0;
    _number_of_frames() = 1;
    _x().setZero(_natoms(),1);
    _y().setZero(_natoms(),1);
    _z().setZero(_natoms(),1);
    _x().block(0,0,mol1._natoms(),1) = mol1._x().col(frame);
    _y().block(0,0,mol1._natoms(),1) = mol1._y().col(frame);
    _z().block(0,0,mol1._natoms(),1) = mol1._z().col(frame);
    _x().block(mol1._natoms(),0,mol2._natoms(),1) = mol2._x().col(frame);
    _y().block(mol1._natoms(),0,mol2._natoms(),1) = mol2._y().col(frame);
    _z().block(mol1._natoms(),0,mol2._natoms(),1) = mol2._z().col(frame);
    return;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sassubset::Mask::
copy_molecule_using_mask(sasmol::SasMol & mol, const boost::dynamic_bitset<> & mask, const size_t frame)
{
    std::vector<int> indices = get_indices_from_mask(mask);
    for (std::vector<int>::size_type i = 0; i<indices.size(); ++i)
    {
        mol._atom_record().push_back(_atom_record()[indices[i]]);
        mol._atom_index().push_back(_atom_index()[indices[i]]);
        mol._atom_name().push_back(_atom_name()[indices[i]]);
        mol._atom_altloc().push_back(_atom_altloc()[indices[i]]);
        mol._atom_resname().push_back(_atom_resname()[indices[i]]);
        mol._atom_chain().push_back(_atom_chain()[indices[i]]);
        mol._atom_resid().push_back(_atom_resid()[indices[i]]);
        mol._atom_icode().push_back(_atom_icode()[indices[i]]);
        mol._atom_occupancy().push_back(_atom_occupancy()[indices[i]]);
        mol._atom_beta().push_back(_atom_beta()[indices[i]]);    
        mol._atom_segname().push_back(_atom_segname()[indices[i]]);    
        mol._atom_element().push_back(_atom_element()[indices[i]]);    
        mol._atom_selement().push_back(_atom_selement()[indices[i]]);    
        mol._atom_charge().push_back(_atom_charge()[indices[i]]);    
    }
    mol._natoms() = indices.size();
    mol._number_of_frames() = 1;
    Eigen::Matrix<float,3,Eigen::Dynamic> coor = get_coor_using_mask(frame, mask);
    mol._x().setZero(mol._natoms(),1);
    mol._y().setZero(mol._natoms(),1);
    mol._z().setZero(mol._natoms(),1);
    mol._x().col(0) = coor.row(0);
    mol._y().col(0) = coor.row(1);
    mol._z().col(0) = coor.row(2);
    return;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
Eigen::Matrix<float,3,Eigen::Dynamic>
sassubset::Mask::
get_coor_using_mask(const size_t frame, const boost::dynamic_bitset<> & mask)
{
    std::vector<int> indices = get_indices_from_mask(mask);
    Eigen::Matrix<float,3,Eigen::Dynamic> coor(3, indices.size());
    for (std::vector<int>::size_type i = 0; i<indices.size(); ++i)
    {
        coor(0,i) = _x().col(frame)[indices[i]];
        coor(1,i) = _y().col(frame)[indices[i]];
        coor(2,i) = _z().col(frame)[indices[i]];
    }
    return coor;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sassubset::Mask::
set_coor_using_mask(sasmol::SasMol & mol, const size_t frame, const boost::dynamic_bitset<> & mask)
{
    std::vector<int> indices = get_indices_from_mask(mask);
    for (std::vector<int>::size_type i = 0; i<indices.size(); ++i)
    {
        dynamic_cast<sasmol::SasMol*>(this)->_x().col(frame)[indices[i]] = mol._x().col(frame)[indices[i]];
        dynamic_cast<sasmol::SasMol*>(this)->_y().col(frame)[indices[i]] = mol._y().col(frame)[indices[i]];
        dynamic_cast<sasmol::SasMol*>(this)->_z().col(frame)[indices[i]] = mol._z().col(frame)[indices[i]];
    }
    return;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sassubset::Mask::
duplicate_molecule(sasmol::SasMol & mol, const int number_of_duplicates, const int frame, std::vector<std::vector<float>>& com_coor, const bool flag_rotate)
{
    const int number_of_frames = 1;
    const int natoms = _natoms();

    std::vector<float> self_com = dynamic_cast<sasmol::SasMol*>(this)->calc_com(const_cast<int&>(frame));

    int count = 1;
    int axis_int;
    std::string axis;
    float theta;
    srand (time(NULL));

    std::string separate_resids = "no";

    mol._x().setZero(_natoms()*number_of_duplicates, 1);
    mol._y().setZero(_natoms()*number_of_duplicates, 1);
    mol._z().setZero(_natoms()*number_of_duplicates, 1);
    
    int i,j;
    std::string this_segment;
    int this_resid=0;
    for (i=0; i<number_of_duplicates; ++i)
    {
        if(i<10) this_segment = std::string("000")+std::to_string(i);
        else if(i<100) this_segment = std::string("00")+std::to_string(i);	
        else if(i<1000) this_segment = std::string("0")+std::to_string(i);
        else this_segment = std::to_string(i);

        if(separate_resids.compare("yes")==0)
        {
            if(i==0) this_resid = 1;
            else ++this_resid;
        }

        for (j=0; j<natoms; ++j)
        {
            mol._atom_record().push_back(_atom_record()[j]);
            mol._atom_index().push_back(count); ++count;
            mol._atom_name().push_back(_atom_name()[j]);
            mol._atom_altloc().push_back(_atom_altloc()[j]);
            mol._atom_resname().push_back(_atom_resname()[j]);
            mol._atom_chain().push_back(_atom_chain()[j]);
            if (separate_resids.compare("yes")==0) mol._atom_resid().push_back(this_resid);
            else mol._atom_resid().push_back(_atom_resid()[j]);
            mol._atom_icode().push_back(_atom_icode()[j]);
            mol._atom_occupancy().push_back(_atom_occupancy()[j]);
            mol._atom_beta().push_back(_atom_beta()[j]);    
            mol._atom_segname().push_back(_atom_segname()[j]);    
            mol._atom_element().push_back(_atom_element()[j]);    
            mol._atom_selement().push_back(_atom_selement()[j]);    
            mol._atom_charge().push_back(_atom_charge()[j]);    
        }

        if (flag_rotate)
        {
            axis_int = rand()%3+1;
            switch(axis_int)
            {
                case 1:
                    axis="x";
                    break;
                case 2:
                    axis="y";
                    break;
                case 3:
                    axis="z";
                    break;
            }
            theta = 360.0 * static_cast<float>(rand()/static_cast<float>(RAND_MAX));
            dynamic_cast<sasmol::SasMol*>(this)->rotate(const_cast<int&>(frame),axis,theta);
        }
        mol._x().col(0).block(_natoms()*i,0,_natoms(),1) = _x().col(frame) + self_com[0];
        mol._y().col(0).block(_natoms()*i,0,_natoms(),1) = _y().col(frame) + self_com[1];
        mol._z().col(0).block(_natoms()*i,0,_natoms(),1) = _z().col(frame) + self_com[2];
        if (flag_rotate) dynamic_cast<sasmol::SasMol*>(this)->rotate(const_cast<int&>(frame),axis,theta);
    }
    mol._natoms() = _natoms()*number_of_duplicates;
    mol._number_of_frames() = 1;
    dynamic_cast<sasmol::SasMol*>(this)->_set_unique_attributes();
}
