      subroutine sub_PullingLD(nres,res_x,res_y,res_z,res_na,
     &     HB_TotNumber,HB_index_I,HB_index_J,TertiaryBond_num,
     &     TertiaryBond_index_I,TertiaryBond_index_J,
     &     TertiaryBond_strength,res_chainid)

      implicit none

      integer nres
      real*8 res_x(nres)
      real*8 res_y(nres)
      real*8 res_z(nres)
  
      character*3 res_na(nres)
      integer HB_TotNumber
      integer HB_index_I(HB_TotNumber)
      integer HB_index_J(HB_TotNumber)
      integer TertiaryBond_num
      integer TertiaryBond_index_I(TertiaryBond_num)
      integer TertiaryBond_index_J(TertiaryBond_num)
      real*8 TertiaryBond_strength(TertiaryBond_num)
      character*1 res_chainid(nres)

      integer nsimu
      parameter (nsimu=1000000)
      real*8 k_bondlength
      parameter (k_bondlength=50.0)
      real*8 k_bondangle_HS
      parameter (k_bondangle_HS=50.0)
      real*8 k_bondangle_L
      parameter (k_bondangle_L=50.0)
      real*8 k_dihedral_HS
      parameter (k_dihedral_HS=50.0)
      real*8 k_dihedral_L
      parameter (k_dihedral_L=0.0)
      real*8 k_disurfide
      parameter (k_disurfide=50.0)
      real*8 k_intradomain
      parameter (k_intradomain=50.0)
      real*8 k_pullingspring
      parameter (k_pullingspring=50.0)
      real*8 cm_frict_const
      parameter (cm_frict_const=10.0)
      real*8 cc_frict_const
      parameter (cc_frict_const=1.0)
      real*8 diffu_const
      parameter (diffu_const=100)
      real*8 dt
      parameter (dt=0.001)
      integer niter
      parameter (niter=100)
      real*8 iter_tol
      parameter (iter_tol=0.0001)
      real*8 min_dist
      parameter (min_dist=1.0)
      real*8 HB_epsilon
      parameter (HB_epsilon=15.0)
      real*8 RP_epsilon,rm_RP,rm_RP_t
      parameter (RP_epsilon=1.0,rm_RP=12.0)
      real*8 Tertiary_weight
      parameter (Tertiary_weight=1.0)
      real*8 HB_LengthCutoff
      parameter (HB_LengthCutoff=2.0)
      integer trj_num
      parameter (trj_num=1)
      real*8 CGHB_angle_cutoff
      parameter (CGHB_angle_cutoff=30.0)
      real*8 ds_distcutoff
      parameter (ds_distcutoff=6.5)
      integer restype2,bin_num
      parameter (restype2=20,bin_num=20)
      real*8 pai
      parameter (pai=3.1415926)

      integer i,j,k
      real*8 x(nres)
      real*8 y(nres)
      real*8 z(nres)
      real*8 vx(nres)
      real*8 vy(nres)
      real*8 vz(nres)
      real*8 rij,dij
      integer num_bondlength
      integer bondlength_pair(nres,2)
      real*8 bondlength_dist(nres)
      integer num_bondangle
      integer bondangle_pair(nres,2)
      real*8 bondangle_dist(nres)
      integer num_dihedral
      integer dihedral_pair(nres,2)
      real*8 dihedral_dist(nres)
      integer num_disurfide
      integer disurfide_pair(nres,2)
      real*8 disurfide_dist(nres)
      integer num_pullingspring
      integer pullingspring_pair(nres)
      real*8 pullingspring_dist(nres)

      integer itime
      real*8 f_bondlength_x(nres)
      real*8 f_bondlength_y(nres)
      real*8 f_bondlength_z(nres)
      real*8 f_bondangle_x(nres)
      real*8 f_bondangle_y(nres)
      real*8 f_bondangle_z(nres)
      real*8 f_dihedral_x(nres)
      real*8 f_dihedral_y(nres)
      real*8 f_dihedral_z(nres)
      real*8 f_disurfide_x(nres)
      real*8 f_disurfide_y(nres)
      real*8 f_disurfide_z(nres)
      real*8 f_HB_x(nres)
      real*8 f_HB_y(nres)
      real*8 f_HB_z(nres)
      real*8 f_SC_x(nres)
      real*8 f_SC_y(nres)
      real*8 f_SC_z(nres)
      real*8 f_RP_x(nres)
      real*8 f_RP_y(nres)
      real*8 f_RP_z(nres)
      real*8 f_rdm_x(nres)
      real*8 f_rdm_y(nres)
      real*8 f_rdm_z(nres)
      real*8 f_pullingspring_x(nres)
      real*8 f_pullingspring_y(nres)
      real*8 f_pullingspring_z(nres)
      integer molecule,neighbor
      real*8 vx_old(nres),vy_old(nres),vz_old(nres)
      real*8 vx_new(nres),vy_new(nres),vz_new(nres)
      real*8 ax(nres,nres),ay(nres,nres),az(nres,nres)
      real*8 bx(nres),by(nres),bz(nres)
      real*8 dir_x,dir_y,dir_z
      real*8 amp_dir
      real*8 scaled_dir_x,scaled_dir_y,scaled_dir_z
      real*8 RMSD
      real*8 pulling_velocity_x1,pulling_velocity_x2
      real*8 pulling_velocity_y1,pulling_velocity_y2
      real*8 pulling_velocity_z1,pulling_velocity_z2
      integer pulling_node1,pulling_node2
      real*8 rm_HB(HB_TotNumber)
      real*8 forceAMP_HB
      real*8 forceAMP_SC
      real*8 forceAMP_RP
      real*8 rm_SC(TertiaryBond_num)
      integer check_flag(nres,nres)
      real*8 RP_R(nres,nres)
      integer itrj
      real*8 Force_endI_trj
      real*8 Force_endJ_trj
      real*8 Force_endI_ave(nsimu)
      real*8 Force_endJ_ave(nsimu)
      real*8 IntEne_ave(nsimu)
      real*8 Hbond_potential(7,36,36)
      real*8 scale
      real*8 local_coor_x_x
      real*8 local_coor_x_y
      real*8 local_coor_x_z
      real*8 local_coor_y_x
      real*8 local_coor_y_y
      real*8 local_coor_y_z
      real*8 local_coor_z_x
      real*8 local_coor_z_y
      real*8 local_coor_z_z
      real*8 mp_x,mp_y,mp_z
      real*8 lx_sc
      real*8 ly_sc
      real*8 lz_sc,l0_sc
      real*8 local_th
      real*8 local_fi
      real*8 local_fi_check
      real*8 dofi_x
      real*8 dofi_y
      real*8 dofi_z
      integer bin_index
      real*8 Hbond_energy
      real*8 doth1,doth2,conv2
      real*8 CGHB_th(HB_TotNumber)
      real*8 CGHB_fi(HB_TotNumber)
      real*8 k_bondangle(nres)
      real*8 k_dihedral(nres)
      real*8 TertiaryBond_strength_2(TertiaryBond_num)
      real*8 FluctRMSD(nres)
      real*8 cm_x,cm_y,cm_z
      integer cm_resnum
      real*8 potential(restype2,restype2,bin_num)
      real*8 ca_potential
      integer dis_ind
      integer i2,j2,index_i,index_j
      real*8 rij_allatom
      integer tempi,tempj
      character*3 res_type(restype2+1)
      integer res_ind(restype2+1)
      real*8 mc1_x,mc1_y,mc1_z
      integer index_mc1
      real*8 mc2_x,mc2_y,mc2_z
      integer index_mc2
	  
      real rand3
      double precision r3

      r3=5.0


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c>>     Multiple Trajectories
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c      print*,'sub'
      do itrj=1,trj_num


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   set pulling velocity
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

         pulling_node1=nres
         
         mc1_x=0
         mc1_y=0
         mc1_z=0
         index_mc1=0

         do i=3309,7340,288
            index_mc1=index_mc1+1
            mc1_x=mc1_x+res_x(i)
            mc1_y=mc1_y+res_y(i)
            mc1_z=mc1_z+res_z(i)
            
         enddo
         mc1_x=mc1_x/real(index_mc1)
         mc1_y=mc1_y/real(index_mc1)
         mc1_z=mc1_z/real(index_mc1)
         
         mc2_x=0
         mc2_y=0
         mc2_z=0
         index_mc2=0
         
         do i=3309,7340
            index_mc2=index_mc2+1
            mc2_x=mc2_x+res_x(i)
            mc2_y=mc2_y+res_y(i)
            mc2_z=mc2_z+res_z(i)
            
         enddo
         mc2_x=mc2_x/real(index_mc2)
         mc2_y=mc2_y/real(index_mc2)
         mc2_z=mc2_z/real(index_mc2)		 
		 
		 
c         open (unit=10,file=
c     &        'tempout.pdb',
c     &        status='unknown',access='append')
c         write(10,2100) 'ATOM  ',1,' CA '
c     &        ,res_na(1), 
c     &        res_chainid(1),1,mc1_x,mc1_y,mc1_z		 
c         write(10,2100) 'ATOM  ',2,' CA '
c     &        ,res_na(2), 
c     &        res_chainid(2),2,mc2_x,mc2_y,mc2_z
c         close(10)
         


         pulling_velocity_x1=(0.0001/dt)*
     &        (mc1_x-mc2_x)/
     &        SQRT((mc2_x-mc1_x)**2+
     &        (mc2_y-mc1_y)**2+
     &        (mc2_z-mc1_z)**2)
	 
         pulling_velocity_y1=(0.0001/dt)*
     &        (mc1_y-mc2_y)/
     &        SQRT((mc2_x-mc1_x)**2+
     &        (mc2_y-mc1_y)**2+
     &        (mc2_z-mc1_z)**2)

         pulling_velocity_z1=(0.0001/dt)*
     &        (mc1_z-mc2_z)/
     &        SQRT((mc2_x-mc1_x)**2+
     &        (mc2_y-mc1_y)**2+
     &        (mc2_z-mc1_z)**2)
	 
c         print*,'pull'


         do i=1,HB_TotNumber
            rm_HB(i)=sqrt((res_x(HB_index_I(i))-
     &           res_x(HB_index_J(i)))**2+
     &           (res_y(HB_index_I(i))-
     &           res_y(HB_index_J(i)))**2+
     &           (res_z(HB_index_I(i))-
     &           res_z(HB_index_J(i)))**2)
         enddo
         
         do i=1,TertiaryBond_num
            rm_SC(i)=sqrt(
     &           (res_x(TertiaryBond_index_I(i))-
     &           res_x(TertiaryBond_index_J(i)))**2+
     &           (res_y(TertiaryBond_index_I(i))-
     &           res_y(TertiaryBond_index_J(i)))**2+
     &           (res_z(TertiaryBond_index_I(i))-
     &           res_z(TertiaryBond_index_J(i)))**2)
         enddo
         
         do i=1,nres
            x(i)=res_x(i)
            y(i)=res_y(i)
            z(i)=res_z(i)
         enddo

c>>  bond length constrain
c         print*,'check1'
         num_bondlength=0
         do i=1,nres-1
            if(res_chainid(i).eq.res_chainid(i+1))then
               rij=sqrt((x(i)-x(i+1))**2
     &              +(y(i)-y(i+1))**2
     &              +(z(i)-z(i+1))**2)
               num_bondlength=num_bondlength+1
               bondlength_dist(num_bondlength)=rij
               bondlength_pair(num_bondlength,1)=i
               bondlength_pair(num_bondlength,2)=i+1
            endif
         enddo
         
         do i=1,nres
            do j=i,nres
               check_flag(i,j)=0
               RP_R(i,j)=0
            enddo
         enddo
         
         do i=1,nres-1
            do j=i+1,nres
               do k=1,TertiaryBond_num
                  if(((TertiaryBond_index_I(k).eq.i)
     &                 .AND.(TertiaryBond_index_J(k).eq.j))
     &                 .OR.((TertiaryBond_index_I(k).eq.j)
     &                 .AND.(TertiaryBond_index_J(k).eq.i)))then
                     check_flag(i,j)=1
                  endif
               enddo
            enddo
         enddo

         do i=1,nres-1
            do j=i+1,nres
               dij=sqrt((x(i)-x(j))**2
     &              +(y(i)-y(j))**2+(z(i)-z(j))**2)            
               if(dij.lt.8.5)then
                  RP_R(i,j)=dij
               else
                  RP_R(i,j)=8.5
               endif
            enddo
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c>>>        main simulation
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         do itime=1,nsimu
            print*,itime

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Force Calculation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c>  calculate the elastic forces between two consecutive c-alpha atom as bond length constrain

            do molecule=1,nres
               f_bondlength_x(molecule)=0
               f_bondlength_y(molecule)=0
               f_bondlength_z(molecule)=0
            enddo
            
            do i=1,num_bondlength
               
               dij=sqrt((x(bondlength_pair(i,1))-
     &              x(bondlength_pair(i,2)))**2+
     &              (y(bondlength_pair(i,1))-
     &              y(bondlength_pair(i,2)))**2+
     &              (z(bondlength_pair(i,1))-
     &              z(bondlength_pair(i,2)))**2)
               
               f_bondlength_x(bondlength_pair(i,1))=
     &              f_bondlength_x(bondlength_pair(i,1))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (x(bondlength_pair(i,2))-
     &              x(bondlength_pair(i,1)))
               
               f_bondlength_y(bondlength_pair(i,1))=
     &              f_bondlength_y(bondlength_pair(i,1))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (y(bondlength_pair(i,2))-
     &              y(bondlength_pair(i,1)))
               
               f_bondlength_z(bondlength_pair(i,1))=
     &              f_bondlength_z(bondlength_pair(i,1))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (z(bondlength_pair(i,2))-
     &              z(bondlength_pair(i,1)))
               
               f_bondlength_x(bondlength_pair(i,2))=
     &              f_bondlength_x(bondlength_pair(i,2))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (x(bondlength_pair(i,1))-
     &              x(bondlength_pair(i,2)))
               
               f_bondlength_y(bondlength_pair(i,2))=
     &              f_bondlength_y(bondlength_pair(i,2))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (y(bondlength_pair(i,1))-
     &              y(bondlength_pair(i,2)))
               
               f_bondlength_z(bondlength_pair(i,2))=
     &              f_bondlength_z(bondlength_pair(i,2))
     &              +k_bondlength*
     &              (1-bondlength_dist(i)/dij)*
     &              (z(bondlength_pair(i,1))-
     &              z(bondlength_pair(i,2)))
               
            enddo
            


c>> calculate the stochastic force of each residue
         
            do molecule=1,nres
               dir_x=2*rand3(r3)-1
               dir_y=2*rand3(r3)-1
               dir_z=2*rand3(r3)-1
               amp_dir=sqrt(dir_x**2+dir_y**2+dir_z**2)
               scaled_dir_x=dir_x/amp_dir
               scaled_dir_y=dir_y/amp_dir
               scaled_dir_z=dir_z/amp_dir
               f_rdm_x(molecule)=5*sqrt(diffu_const)*  ! original value is 2
     &              cm_frict_const*scaled_dir_x
               f_rdm_y(molecule)=5*sqrt(diffu_const)*
     &              cm_frict_const*scaled_dir_y
               f_rdm_z(molecule)=5*sqrt(diffu_const)*
     &              cm_frict_const*scaled_dir_z
            enddo

c>  calculate the inter-domain Hydrogen bond as LJ potential; but intradomain as harmonic

            do molecule=1,nres
               f_HB_x(molecule)=0
               f_HB_y(molecule)=0
               f_HB_z(molecule)=0
            enddo

            do i=1,HB_TotNumber

               dij=sqrt((x(HB_index_I(i))-
     &              x(HB_index_J(i)))**2+
     &              (y(HB_index_I(i))-
     &              y(HB_index_J(i)))**2+
     &              (z(HB_index_I(i))-
     &              z(HB_index_J(i)))**2)
               

               f_HB_x(HB_index_I(i))=
     &              f_HB_x(HB_index_I(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (x(HB_index_J(i))-
     &              x(HB_index_I(i)))
               f_HB_y(HB_index_I(i))=
     &              f_HB_y(HB_index_I(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (y(HB_index_J(i))-
     &              y(HB_index_I(i)))
               f_HB_z(HB_index_I(i))=
     &              f_HB_z(HB_index_I(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (z(HB_index_J(i))-
     &              z(HB_index_I(i)))
               
               f_HB_x(HB_index_J(i))=
     &              f_HB_x(HB_index_J(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (x(HB_index_I(i))-
     &              x(HB_index_J(i)))
               f_HB_y(HB_index_J(i))=
     &              f_HB_y(HB_index_J(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (y(HB_index_I(i))-
     &              y(HB_index_J(i)))
               f_HB_z(HB_index_J(i))=
     &              f_HB_z(HB_index_J(i))
     &              +k_intradomain*
     &              (1-rm_HB(i)/dij)*
     &              (z(HB_index_I(i))-
     &              z(HB_index_J(i)))
               

            enddo
   
c>  calculate the interdomain Side-Chain Packing as LJ potential; but intradomain as Harmonic

            do molecule=1,nres
               f_SC_x(molecule)=0
               f_SC_y(molecule)=0
               f_SC_z(molecule)=0
            enddo
            
            do i=1,TertiaryBond_num
               
               dij=sqrt((x(TertiaryBond_index_I(i))-
     &              x(TertiaryBond_index_J(i)))**2+
     &              (y(TertiaryBond_index_I(i))-
     &              y(TertiaryBond_index_J(i)))**2+
     &              (z(TertiaryBond_index_I(i))-
     &              z(TertiaryBond_index_J(i)))**2)


               f_SC_x(TertiaryBond_index_I(i))=
     &              f_SC_x(TertiaryBond_index_I(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*
     &              (x(TertiaryBond_index_J(i))-
     &              x(TertiaryBond_index_I(i)))
               f_SC_y(TertiaryBond_index_I(i))=
     &              f_SC_y(TertiaryBond_index_I(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*
     &              (y(TertiaryBond_index_J(i))-
     &              y(TertiaryBond_index_I(i)))
               f_SC_z(TertiaryBond_index_I(i))=
     &              f_SC_z(TertiaryBond_index_I(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*
     &              (z(TertiaryBond_index_J(i))-
     &              z(TertiaryBond_index_I(i)))
               
               f_SC_x(TertiaryBond_index_J(i))=
     &              f_SC_x(TertiaryBond_index_J(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*            
     &              (x(TertiaryBond_index_I(i))-
     &              x(TertiaryBond_index_J(i)))
               f_SC_y(TertiaryBond_index_J(i))=
     &              f_SC_y(TertiaryBond_index_J(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*
     &              (y(TertiaryBond_index_I(i))-
     &              y(TertiaryBond_index_J(i)))
               f_SC_z(TertiaryBond_index_J(i))=
     &              f_SC_z(TertiaryBond_index_J(i))
     &              +k_intradomain*
     &              (1-rm_SC(i)/dij)*
     &              (z(TertiaryBond_index_I(i))-
     &              z(TertiaryBond_index_J(i)))
               
            enddo

 
c>  calculate the repulsive interaction between residues

            do molecule=1,nres
               f_RP_x(molecule)=0
               f_RP_y(molecule)=0
               f_RP_z(molecule)=0
            enddo
            
c            do i=1,nres-1
c               do j=i+1,nres

c                  if((res_chainid(i).eq.'R').and.
c     &                 (res_chainid(j).eq.'R'))then
                     
c                     if((abs(i-j).ge.2).and.
c     &                    (check_flag(i,j).eq.0))then

c                        dij=sqrt((x(i)-x(j))**2
c     &                       +(y(i)-y(j))**2
c     &                       +(z(i)-z(j))**2)
                     
c                        if(dij.le.20.0)then

c                           rm_RP_t=RP_R(i,j)

c                           forceAMP_RP=-6*RP_epsilon*(1/dij)
c     &                          *((rm_RP_t/dij)**6)*
c     &                          (2.0+2.0*((rm_RP_t/dij)**6))
                           
c                           f_RP_x(i)=f_RP_x(i)+
c     &                          forceAMP_RP*(x(j)-x(i))/dij
c                           f_RP_y(i)=f_RP_y(i)+
c     &                          forceAMP_RP*(y(j)-y(i))/dij
c                           f_RP_z(i)=f_RP_z(i)+
c     &                          forceAMP_RP*(z(j)-z(i))/dij
c                           
c                           f_RP_x(j)=f_RP_x(j)
c     &                          +forceAMP_RP*(x(i)-x(j))/dij
c                           f_RP_y(j)=f_RP_y(j)
c     &                          +forceAMP_RP*(y(i)-y(j))/dij
c                           f_RP_z(j)=f_RP_z(j)
c     &                          +forceAMP_RP*(z(i)-z(j))/dij
c                           
c                        endif
c                        
c                     endif

c                  endif
                       
                  
c               enddo
c            enddo


            do i=1,nres-1
               do j=i+1,nres

                  if(((res_chainid(i).eq.'R').and.
     &                 (res_chainid(j).ne.'R')).or.
     &                 ((res_chainid(i).ne.'R').and.
     &                 (res_chainid(j).eq.'R')))then

                     dij=sqrt((x(i)-x(j))**2
     &                    +(y(i)-y(j))**2
     &                    +(z(i)-z(j))**2)
                     
                     if(dij.le.8.5)then

                        rm_RP_t=RP_R(i,j)
                           
                        forceAMP_RP=-6*RP_epsilon*(1/dij)
     &                       *((rm_RP_t/dij)**6)*
     &                       (2.0+2.0*((rm_RP_t/dij)**6))
                        
                        f_RP_x(i)=f_RP_x(i)+
     &                       forceAMP_RP*(x(j)-x(i))/dij
                        f_RP_y(i)=f_RP_y(i)+
     &                       forceAMP_RP*(y(j)-y(i))/dij
                        f_RP_z(i)=f_RP_z(i)+
     &                       forceAMP_RP*(z(j)-z(i))/dij
                        
                        f_RP_x(j)=f_RP_x(j)
     &                       +forceAMP_RP*(x(i)-x(j))/dij
                        f_RP_y(j)=f_RP_y(j)
     &                       +forceAMP_RP*(y(i)-y(j))/dij
                        f_RP_z(j)=f_RP_z(j)
     &                       +forceAMP_RP*(z(i)-z(j))/dij
                           
                     endif

                  endif
                       
                  
               enddo
            enddo
           
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>> calculate the velocity of each cell
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            do molecule=1,nres
               vx_old(molecule)=0
               vy_old(molecule)=0
               vz_old(molecule)=0
               do neighbor=1,nres
                  ax(molecule,neighbor)=0
                  ay(molecule,neighbor)=0
                  az(molecule,neighbor)=0
               enddo
               bx(molecule)=0
               by(molecule)=0
               bz(molecule)=0
            enddo
            
            do molecule=1,nres
               ax(molecule,molecule)=
     &              ax(molecule,molecule)+cm_frict_const
               ay(molecule,molecule)=
     &              ay(molecule,molecule)+cm_frict_const
               az(molecule,molecule)=
     &              az(molecule,molecule)+cm_frict_const
               do neighbor=1,nres
                  if(molecule.ne.neighbor)then
                     dij=sqrt((x(molecule)-x(neighbor))**2+
     &                    (y(molecule)-y(neighbor))**2+
     &                    (z(molecule)-z(neighbor))**2)
                     if(dij.lt.min_dist)then
                        ax(molecule,molecule)=
     &                       ax(molecule,molecule)+
     &                       cc_frict_const
                        ax(molecule,neighbor)=
     &                       ax(molecule,neighbor)+
     &                       cc_frict_const
                        ay(molecule,molecule)=
     &                       ay(molecule,molecule)+
     &                       cc_frict_const
                        ay(molecule,neighbor)=
     &                       ay(molecule,neighbor)+
     &                       cc_frict_const
                        az(molecule,molecule)=
     &                       az(molecule,molecule)+
     &                       cc_frict_const
                        az(molecule,neighbor)=
     &                       az(molecule,neighbor)+
     &                       cc_frict_const
                     endif
                  endif
               enddo
            enddo
            
            do molecule=1,nres
               bx(molecule)=f_bondlength_x(molecule)
     &              +f_rdm_x(molecule)
     &              +f_HB_x(molecule)
     &              +f_SC_x(molecule)
     &              +f_RP_x(molecule)
               by(molecule)=f_bondlength_y(molecule)
     &              +f_rdm_y(molecule)
     &              +f_HB_y(molecule)
     &              +f_SC_y(molecule)
     &              +f_RP_y(molecule)
               bz(molecule)=f_bondlength_z(molecule)
     &              +f_rdm_z(molecule)
     &              +f_HB_z(molecule)
     &              +f_SC_z(molecule)
     &              +f_RP_z(molecule)
            enddo
            
            do i=1,niter
               
               do molecule=1,nres
                  vx_new(molecule)=bx(molecule)
                  vy_new(molecule)=by(molecule)
                  vz_new(molecule)=bz(molecule)
                  do neighbor=1,nres
                     if(molecule.ne.neighbor)then
                        vx_new(molecule)=vx_new(molecule)+
     &                       ax(molecule,neighbor)*vx_old(neighbor)
                        vy_new(molecule)=vy_new(molecule)+
     &                       ay(molecule,neighbor)*vy_old(neighbor)
                        vz_new(molecule)=vz_new(molecule)+
     &                       az(molecule,neighbor)*vz_old(neighbor)
                     endif
                  enddo
                  vx_new(molecule)=vx_new(molecule)
     &                 /ax(molecule,molecule)
                  vy_new(molecule)=vy_new(molecule)
     &                 /ay(molecule,molecule)
                  vz_new(molecule)=vz_new(molecule)
     &                 /az(molecule,molecule)
               enddo
               
               RMSD=0
               do molecule=1,nres
                  RMSD=RMSD+(vx_new(molecule)-vx_old(molecule))**2+
     &                 (vy_new(molecule)-vy_old(molecule))**2+
     &                 (vz_new(molecule)-vz_old(molecule))**2
               enddo
               RMSD=sqrt(RMSD/real(nres))
               
               if(RMSD.lt.iter_tol)then
                  goto 100
               endif
               do molecule=1,nres
                  vx_old(molecule)=vx_new(molecule)
                  vy_old(molecule)=vy_new(molecule)
                  vz_old(molecule)=vz_new(molecule)
               enddo
               
            enddo
            
 100        continue
            
            do molecule=1,nres
               vx(molecule)=vx_new(molecule)
               vy(molecule)=vy_new(molecule)
               vz(molecule)=vz_new(molecule)
            enddo



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>  update the position of each cell
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            
            do j=1,nres
               if(j.le.572)then
                  x(j)=x(j)
     &                 +pulling_velocity_x1*dt
                  y(j)=y(j)
     &                 +pulling_velocity_y1*dt
                  z(j)=z(j)
     &                 +pulling_velocity_z1*dt			   
               else
                  x(j)=x(j)+vx(j)*dt
                  y(j)=y(j)+vy(j)*dt
                  z(j)=z(j)+vz(j)*dt   
               endif 
            enddo
                        

cccccccccccccccccccccccccccccccccccccccccccccc
c>>>   data output
cccccccccccccccccccccccccccccccccccccccccccccc


            if(MOD(itime,100).eq.1)then
               open (unit=10,file=
     &              'PullingPP7pilus_trj_k50v0001.pdb',
     &              status='unknown',access='append')
               do j=1,nres
                  if(res_chainid(j).ne.'R')then
                     write(10,2100) 'ATOM  ',j,' CA '
     &                    ,res_na(j), 
     &                    res_chainid(j),j,x(j),y(j),z(j)
                  endif
               enddo
               write(10,2102) 'TER'
               do j=1,nres
                  if(res_chainid(j).eq.'R')then
                     write(10,2100) 'ATOM  ',j,' CA '
     &                    ,res_na(j), 
     &                    res_chainid(j),j,x(j),y(j),z(j)
                  endif
               enddo
               write(10,2102) 'TER'
               write(10,2102) 'END'
               close(10)
               
            endif

 2100       format(A6,I5,1x,A4,1x,A3,1x,A1,I4,4x,3F8.3)
 2102       format(A3)      
       


         enddo

      enddo


      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
