      program main
      implicit none

      integer nres_C
      parameter (nres_C=3308)
      integer nres_R
      parameter (nres_R=4032)
      integer nres
      parameter (nres=nres_C+nres_R)


      character*3 res_na(nres)
      character*1 chainid(nres)
      real*8 x(nres)
      real*8 y(nres)
      real*8 z(nres)

      integer i,j,index
      integer HB_TotNumber
      integer HB_index_I(nres*20)
      integer HB_index_J(nres*20)
      integer TertiaryBond_num
      integer TertiaryBond_index_I(nres*100)
      integer TertiaryBond_index_J(nres*100)
      real*8 TertiaryBond_strength(nres*100)
      real*8 dij
      integer i2,j2,T_flag
      integer seq_bg,seq_ed


      open(unit=10,file='PP7_pilus_CA.pdb',
     &     status='old')

      index=0
      do j=1,nres_C
         index=index+1
         read(10,2000) res_na(index),
     &        x(index),y(index),z(index)
         chainid(index)='C'
      enddo
c      read(10,*)
       do j=1,nres_R
         index=index+1
         read(10,2000) res_na(index),
     &        x(index),y(index),z(index)
         chainid(index)='R'	 
      enddo
      close(10)
c      print*,index
 2000 format(17x,A3,10x,3F8.3)


      HB_TotNumber=0
      do i=1,nres
         HB_index_I(i)=0
         HB_index_J(i)=0
      enddo
      
      do i=1,nres-4
         do j=i+4,nres
            if(((chainid(i).eq.'C').and.(chainid(j).eq.'C')).or.
     &           ((chainid(i).eq.'R').and.(chainid(j).eq.'R')))then
               dij=sqrt((x(i)-x(j))**2
     &              +(y(i)-y(j))**2+(z(i)-z(j))**2)
               if(dij.lt.7.5)then
                  HB_TotNumber=HB_TotNumber+1    
                  HB_index_I(HB_TotNumber)=i
                  HB_index_J(HB_TotNumber)=j
c               print*,HB_TotNumber,
c     &              HB_index_I(HB_TotNumber),
c     &              HB_index_J(HB_TotNumber)
               endif
            endif
         enddo
      enddo


      TertiaryBond_num=0
      do i=1,nres
         TertiaryBond_index_I(i)=0
         TertiaryBond_index_J(i)=0
         TertiaryBond_strength(i)=0
      enddo
    
      do i=1,nres-2
         do j=i+2,nres
            T_flag=0            
            dij=sqrt((x(i)-x(j))**2
     &           +(y(i)-y(j))**2
     &           +(z(i)-z(j))**2)
            if(dij.lt.8.5)then
               if(((chainid(i).eq.'C').and.(chainid(j).eq.'C')).or.
     &              ((chainid(i).eq.'R').and.(chainid(j).eq.'R')))then
                  T_flag=1
               endif
            endif
           
            if(T_flag.eq.1)then
               TertiaryBond_num=TertiaryBond_num+1
               TertiaryBond_index_I(TertiaryBond_num)=i
               TertiaryBond_index_J(TertiaryBond_num)=j
               TertiaryBond_strength(TertiaryBond_num)=1
c               print*,TertiaryBond_num,
c     &                TertiaryBond_index_I(TertiaryBond_num),
c     &                TertiaryBond_index_J(TertiaryBond_num)
            endif
         enddo
      enddo

c      print*,'cell'
      call sub_PullingLD(nres,x,y,z,res_na,
     &     HB_TotNumber,HB_index_I,HB_index_J,TertiaryBond_num,
     &     TertiaryBond_index_I,TertiaryBond_index_J,
     &     TertiaryBond_strength,chainid)

      stop
      end
