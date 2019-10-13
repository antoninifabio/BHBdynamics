      SUBROUTINE read()
      USE commonV
      USE modelPar
      USE inputPar
      implicit none
      integer j,i,i2,h,k,l,m,n,kick
      real*8 fbh0,mr,Mtot,psi0,rh_i,Mcl_i,pi,rho,vesc  
      logical :: file_exists,mbh_exists
      pi=2.d0*ASIN(1.d0)    

!     generate simple cluster model if cluster.txt is not found or input=0
     
      INQUIRE(FILE="cluster.txt", EXIST=file_exists)
      if(file_exists.eqv..False.)then
         print*,'no cluster.txt file found'
         input=0
      end if

      if(input.eq.0)then              
         print*,'SimplE cluster model will be generated'
         write(*,*)'Cluster mass (M_Sun), and half-mass radius (pc):'
         read(*,*)Mcl_i,rh_i
         write(*,*)'initial mass fraction in BHs:'
         read(*,*)fbh0
         call simple(Nmax,rh_i,Mcl_i,fbh0) !IMF=1 (0) read (build) the m_max(M)     
         input=1
      end if 
      
!     read/build data      
      inp: select case(input)            
      case(1)                   !cluster ev. + m_max(t) in cluster.txt
!     store cluster properties 
         open (unit=88,file="cluster.txt",action="READ",status="OLD") !file with: time, M_bh, M_cl, r_h, m_max         
         i=0
         do            
            i=i+1
            read(88,*,end=11)tabt(i),tabmbh(i)
     &           ,tabm(i),tabr(i),tabmax(i) !read cluster evolution                        
         end do         
 11      cont=i
         
         if(IMF.eq.1)CALL FMMsplineDouble(tabt,tabmax,c10,c11,c12,cont)         
      case(2)                   !cluster ev. + BH masses in cluster.txt and mbh.txt
!     store cluster properties 
         open (unit=88,file="cluster.txt",action="READ",status="OLD") !file with: time, M_bh, M_cl, r_h
         i=0
         do  
            i=i+1
            read(88,*,end=12)tabt(i),tabmbh(i),tabm(i),tabr(i) !read cluster evolution
         end do
 12      cont=i
!     build M(<m)
         if(IMF.eq.1)then       !skip for a power law IMF
            INQUIRE(FILE="mbh.txt", EXIST=mbh_exists)
            if(mbh_exists.eqv..False.)then
               print*,'mbh.txt not found'
               stop
            end if
            open (unit=78,file="mbh.txt",action="READ",status="OLD") !file with masses of BHs after SN kicks          
            j=0
            mr=0.d0
            do !while (mr.lt.tabmbh(1)) !read BH masses
               j=j+1
               read(78,*,end=13)w(j)
               mr=mr+w(j)            
            end do
 13         call SORTRX(j,w,INDEX) !sort m()
            mr=tabmbh(1)/mr  
            
            Mtot=0.d0
            do i=1,j            
               tabmax(i)=w(INDEX(i))            
               Mtot=tabmax(i)+Mtot
               tabmbh_p(i)=Mtot*mr            
            end do
            cont2=i         
            CALL FMMsplineDouble(tabmbh_p,tabmax,c10,c11,c12,cont2) !make M(<m)
         end if
         
      case DEFAULT
         write(*,*)"value of input not valid"
         stop
      end select inp
      
      if(tabt(cont-1).lt.tf/1.d9.and.tabmbh(cont-1).gt.1)then
         print*,tabt(cont-1),tf,
     &        "warning: cluster simulation not finished"         
      end if
       
!     spline the cluster model
      CALL FMMsplineDouble(tabt,tabmbh,c1,c2,c3,cont)
      CALL FMMsplineDouble(tabt,tabm,c4,c5,c6,cont)
      CALL FMMsplineDouble(tabt,tabr,c7,c8,c9,cont)

      if(tabmbh(1)/tabm(1).gt.0.1)then
         print*,"reduce mass fraction in BHs below 10 per cent"
         stop
      end if

      rho=3.d0*tabm(1)/(8.d0*pi*tabr(1)**3)
      vesc=50.d0*(tabm(1)/1.d5)**(1./3.)*(rho/1.d5)**(1./6.)
      if(vesc.gt.300.and.rho.gt.1.e5)then
         print*,
     &    "vesc>300km/s; hierarchical growth expected" !hierarchical growth expected (Antonini et al.'18)
         stop
      end if

!     compute time when balanced evolution starts (AG'19)
      call TBalanced(tabmbh(1),tabm(1),tabr(1),tcc)      
      
      RETURN
      END Subroutine
