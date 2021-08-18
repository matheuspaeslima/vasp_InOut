It is an oriented object Fortran code for handling VASP output files. 
The handling subroutines are stored in the file m_contcar.f
The user should modify the file "main.f" for a specific use. 
For example, for conversion of cartesian to direct coordinates, the user can modify de "main.f" file as follows:

=====================================
c  .: reading initial geometry    
          fname="POSCAR"
          call read_contcar(data1,fname)
c .: cartesian to direct coordinate transformation:
          call direct_to_cartesian(data1)
c.: writing the final POSCAR
         call write_poscar(data1)
=======================================

The available operations follow above: 


c    0) type(contcar)                (06/04/2017)
c    1) read_contcar(var,fname)      (06/04/2017)
c    2) CONTCAR_to_xyz(var)          (06/04/2017)                           
c    3) write_contcar(var,fname)     (07/04/2017)                     
c    4) direct_to_Cartesian(var)     (07/04/2017)                                    
c    5) Cartesian_to_direct(var)     (07/04/2017)
c    6) translate                    (07/04/2017)                
c    7) put_selective(var)           (07/04/2017)
c    8) put_in_box(var)              (07/04/2017)
c    9) erease_atom(var,i)           (10/04/2017)
c   10) supercell(var,nx,ny,nz)      (10/04/2017)
c   11) locpot skills                (19/06/2017)
c   12) supercellM                   (25/07/2018)
c   13) RotateCell                   (25/07/2018)
c   14) Rotate                       (21/03/2019)
c   15) add_atom                     (01/08/2020)
c   16) purge_contcar                (11/08/2020)
c   17) write_contcarF               (11/08/2020)
c   18) CONTCAR_to_xyzF              (11/08/2020) 
