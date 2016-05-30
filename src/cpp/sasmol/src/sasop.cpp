#include <sasop.h>
#include <sasmol.h>

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasop::Move::
translate(int &frame, std::vector<float> &value){

    _x().col(frame) += value[0] ;
    _y().col(frame) += value[1] ;
    _z().col(frame) += value[2] ;
    
    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasop::Move::
center(int &frame)
{
     std::vector<float> com = dynamic_cast<sasmol::SasMol*>(this)->calc_com(frame) ;

     _x().col(frame) -= com[0] ;
     _y().col(frame) -= com[1] ;
     _z().col(frame) -= com[2] ;
    
     return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasop::Move::
move_to(int &frame, std::vector<float> &value)
{
    std::vector<float> com = dynamic_cast<sasmol::SasMol*>(this)->calc_com(frame) ;

    _x().col(frame) -= (com[0] + value[0]) ;
    _y().col(frame) -= (com[1] + value[1]) ;
    _z().col(frame) -= (com[2] + value[2]) ;
    
    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasop::Move::
align(sasmol::SasMol &mol, int &frame){

    /* NEED TO PASS SAME SUBSET ARGUMENTS FOR CORRECT COM ETC.
    
    NOT DEBUGGED OR CHECKED AT ALL

        Alignment of one object on top of another
        "self" is aligned onto "other" using the basis
        of molecule 2 to align onto the basis of molecule 1
        and the transformation is then done to all the atoms of
        molecule 2
    */

    //Eigen::Matrix3f u = find_u(mol,frame) ;
    Eigen::Matrix3f u = mol.find_u(*dynamic_cast<sasmol::SasMol*>(this),frame) ;

    std::cout << "in align: u = \n" << u << std::endl ;

    std::vector<float> com1 = dynamic_cast<sasmol::SasMol*>(this)->calc_com(frame) ;
    std::vector<float> com2 = mol.calc_com(frame) ;

    // this com needs to be correct !
    
    Eigen::Matrix<float,3,Eigen::Dynamic> coor(3,_natoms()) ; //  (rows,cols) 

    coor.row(0) = _x().col(frame) - com1[0] ; 
    coor.row(1) = _y().col(frame) - com1[1] ;
    coor.row(2) = _z().col(frame) - com1[2] ;                            
    
    Eigen::Matrix<float,3,Eigen::Dynamic> ncoor(3,_natoms()) ; 

    ncoor = (u * coor) ;

    _x().col(frame) = ncoor.row(0).array() + com2[0] ;
    _y().col(frame) = ncoor.row(1).array() + com2[1] ;
    _z().col(frame) = ncoor.row(2).array() + com2[2] ;
    
    return ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
void
sasop::Move::
rotate(int &frame, std::string &axis, float &theta)
{
    /*
        Simple rotation about the x, y, or z axis.
        Note that calcuations are in radians
    */

    float cs = cos(theta), si = sin(theta) ;
    Eigen::Matrix3f mat ;

    if(axis == "x")
    {
        mat << 1.0,0.0,0.0,0.0,cs,-si,0.0,si,cs ;
    }
    else if(axis == "y")
    {
        mat << cs,0.0,si,0.0,1.0,0.0,-si,0.0,cs ;
    }
    else if(axis == "z")
    {
        mat << cs,-si,0.0,si,cs,0.0,0.0,0.0,1.0 ;
    }
    else
    {
        mat.setIdentity() ; std::cout << "need to supply x,y, or z axis\n" ;
    }

    //std::cout << "mat = \n" << mat << std::endl ;

    Eigen::Matrix<float,3,Eigen::Dynamic> coor(3,_natoms()) ; //  (rows,cols) 

    //    coor.rows() = 3 ; coor.cols() = 6730 ;

    // _x().rows() = 6730 --> natoms of coordinates row-major vector // 0
     // _x().cols() = 1   --> frame                        // 1
                                                // 2
                                                // .
                                                // .
                                                // .
                                                // natoms - 1
    // (rows,cols)
    //
    // ( 3 x 1 ) = ( 3 x 3 ) * (3 x 1) 
    //
    //    | x'|     = R * | x |
    //    | y'|           | y | 
    //    | z'|        | z |
    // 

    std::vector<float> com = dynamic_cast<sasmol::SasMol*>(this)->calc_com(frame) ;
    //std::vector<float> com = { 0.0, 0.0, 0.0 } ;

    //  all from a single column to first row  // 0 1 2 . . . natoms - 1 

    coor.row(0) = _x().col(frame) - com[0] ; 
    coor.row(1) = _y().col(frame) - com[1] ;
    coor.row(2) = _z().col(frame) - com[2] ;                            
    
    Eigen::Matrix<float,3,Eigen::Dynamic> ncoor(3,_natoms()) ; 
    ncoor = mat * coor ;
        // (3 x 6730) = (3 x 3) * (3 x 6730) 

    _x().col(frame) = ncoor.row(0).array() + com[0] ;
    _y().col(frame) = ncoor.row(1).array() + com[1] ;
    _z().col(frame) = ncoor.row(2).array() + com[2] ;

/*
eigen (rows,cols)  ... column major is default (which is fortran style)
    numpy is row major (which is C style)
*/


    return ;
}


/*
to do:
    def align(self,frame,coor_sub_2,com_sub_2,coor_sub_1,com_sub_1):

        '''

        '''
        self.masscheck(frame)
        self.calccom(frame)
    
        u = sasmath.find_u(coor_sub_1, coor_sub_2)

        tao = numpy.transpose(self.coor()[frame] - com_sub_2)

        error,nat2 = sasmath.matrix_multiply(u,tao)

        ncoor = numpy.transpose(nat2) + com_sub_1

        self._coor[frame,:] = ncoor
 
        return

*/


