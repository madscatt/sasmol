#include <sasmath.h>
#include <sasmol.h>


/*             u = sasmath.find_u(coor_sub_1, coor_sub_2) */

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
Eigen::Matrix3f
sasmath::Math::
find_u(sasmol::SasMol &mol, int &frame)
{

    Eigen::Matrix3f r ;

    float rxx = (mol._x().col(frame)*_x().col(frame)).sum() ; 
    float rxy = (mol._x().col(frame)*_y().col(frame)).sum() ;
    float rxz = (mol._x().col(frame)*_z().col(frame)).sum() ;
    
    float ryx = (mol._y().col(frame)*_x().col(frame)).sum() ;
    float ryy = (mol._y().col(frame)*_y().col(frame)).sum() ;
    float ryz = (mol._y().col(frame)*_z().col(frame)).sum() ;

    float rzx = (mol._z().col(frame)*_x().col(frame)).sum() ;
    float rzy = (mol._z().col(frame)*_y().col(frame)).sum() ;
    float rzz = (mol._z().col(frame)*_z().col(frame)).sum() ;

    r << rxx,rxy,rxz,
         ryx,ryy,ryz,
         rzx,rzy,rzz;

    Eigen::Matrix3f rtr = r.transpose() * r ;
    
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(rtr);

    if (eigensolver.info() != Eigen::Success)
    {
          std::cout << "find_u calculation failed" << std::endl ;
    }

    Eigen::Matrix<float,3,1> uk ;  // eigenvalues
    Eigen::Matrix3f ak ;          // eigenvectors

    uk = eigensolver.eigenvalues() ; 
    ak = eigensolver.eigenvectors() ;

    Eigen::Matrix3f akt = ak.transpose() ;
    Eigen::Matrix3f new_ak  ;

    new_ak.row(0) = akt.row(2) ; //sort eigenvectors
    new_ak.row(1) = akt.row(1) ;
    new_ak.row(2) = akt.row(0) ;

    // python ak0 --> ak[2] ; ak1 --> ak[1] ; ak2 --> ak[0]

    //Eigen::Matrix<float,3,1> ak2 = akt.row(2).cross(akt.row(1)) ;

    Eigen::MatrixXf rak0 = (r * (akt.row(2).transpose())) ;
    Eigen::MatrixXf rak1 = (r * (akt.row(1).transpose())) ;

    Eigen::Matrix<float,3,1> urak0 ; 
    if(uk[2] == 0.0)
    { 
        urak0 = 1E15 * rak0 ;
    } 
    else
    {    
        urak0 = (1.0/sqrt(fabs(uk[2])))*rak0 ;
    }
    
    Eigen::Matrix<float,3,1> urak1 ; 
    if(uk[1] == 0.0)
    { 
        urak1 = 1E15 * rak1 ;
    } 
    else
    {    
        urak1 = (1.0/sqrt(fabs(uk[1])))*rak1 ;
    }

    Eigen::Matrix<float,3,1> urak2 = urak0.col(0).cross(urak1.col(0)) ;
    Eigen::Matrix3f bk ;

    bk << urak0,urak1,urak2;

    Eigen::Matrix3f u ;

    u =  (bk * new_ak).transpose() ;

    return u ;
/*
u = array([[-0.94513198,  0.31068658,  0.100992  ],
       [-0.3006246 , -0.70612572, -0.64110165],
       [-0.12786863, -0.63628635,  0.76078203]])
check:

u = 
-0.945133  0.310686  0.100992
-0.300624 -0.706126 -0.641102
-0.127869 -0.636287  0.760782 
*/
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sasmath::
signed_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3)
{

    /*
    This method calcultes the sign of the angle which is used in the calculation of a dihedral angle.
    As such, this method is not validated for other uses.  It will fail if the basis atoms for the
    dihedral (atom 2 and atom 3) overlap.
    */

    float sa = 180.0 ;
    float angle ; 
    float pi = acos(-1.0) ;

    float ada = coor1.dot(coor1) ;    
    float bdb = coor2.dot(coor2) ;    

    if(ada*bdb <= 0.0)
    {
        std::cout << "; " ;
        return sa ;
    }
    else
    {
        float adb = coor1.dot(coor2) ;    
        float argument = adb/sqrt(ada*bdb) ;
        angle = (180.0/pi) * acos(argument) ;
    }

    Eigen::Matrix<float,3,1> cp ; cp = coor1.cross(coor2) ;
    float dp = coor3.dot(cp) ;
    
    float sign = 0.0 ;
    if(dp < 0.0) 
    {
        sign = -1.0 ;
    }
    else if(dp > 0.0)
    {
        sign = 1.0 ;
    }

    sa = sign*angle ;

    return sa ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sasmath::
dihedral_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3, Eigen::Matrix<float,3,1> &coor4)
{

    float da ;

    Eigen::Matrix<float,3,1> r1 ; r1 = coor1 - coor2 ;
    Eigen::Matrix<float,3,1> r2 ; r2 = coor3 - coor2 ;
    Eigen::Matrix<float,3,1> r3 ; r3 = -1.0*r2 ;
    Eigen::Matrix<float,3,1> r4 ; r4 = coor4 - coor3 ;
    
    Eigen::Matrix<float,3,1> n1 ; n1 = r1.cross(r2) ;
    Eigen::Matrix<float,3,1> n2 ; n2 = r3.cross(r4) ;

    da = signed_angle(n1,n2,r2) ;

    return da ;
}

///
/// @par Detailed description 
/// ... 
/// @param [in, out] (param1) ...
/// @return ...
/// @note ...
float
sasmath::
calc_angle(Eigen::Matrix<float,3,1> &coor1, Eigen::Matrix<float,3,1> &coor2, Eigen::Matrix<float,3,1> &coor3)
{
    /*
    Method to calculate the angle between three atoms
    */

    float angle ;

    Eigen::Matrix<float,3,1> u ; u = coor1 - coor2 ;
    Eigen::Matrix<float,3,1> v ; v = coor3 - coor2 ;

    float c = (u.dot(v))/(u.norm()*v.norm()) ;
    angle = acos(c) ;

    return angle ;
}


