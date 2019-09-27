!     BHBdynamics computes the evolution of a cluster model, 
!     the number of BH mergers it produces vs time, and 
!     their (cumulative) eccentricity distributions. 
!     See Antonini and Gieles (2019) for reference      
      PROGRAM BHBdynamics
      USE modelPar
      USE inputPar
      implicit none
      real*8, external:: fun,funej2,funej     
      real*8 abserr,dec,dt,eccen,epsabs,epsrel,Fej1,Fej2,tfin
      real*8 Fej,Fin,Ntot,step,t,td,td0,X0,X1,Gamma,N0,N1,tf0
      integer j,i,i2,h,k,l,m,n,cont,cont2,kick   
      common/ecc/eccen,X0,X1
      
!     parameters of integration
      integer limit,lenw,key,ier
      parameter ( limit = 1000 )	
      parameter ( lenw = 4000 )
      parameter ( key = 6)
      parameter ( epsabs =0.0d0)
      parameter ( epsrel =1.0d-6) !integral accuracy
      integer iwork(limit)
      double precision work(lenw)
      integer last
      integer neval 
!     parameters of integration
       
!     model parameters   
      csi=0.1d0                 !relaxation coeff.
      a1=1.91                   !coefficient relating f_bh to t_relax (see Antonini+Gieles)
      m_mean=0.638d0            !mean mass
      beta=2.55d-3              !BH loss rate      
      dE=0.2d0                  !fractional energy change per interaction
!     model parameters

      open (unit=20,file="input.ini",action="READ",status="OLD")
      read(20,*)
      read(20,*)input,IMF,a,tf,td,freq
      tf=tf*1.d9                !cluster formation time
      td=td*1.d9                !look-back time of detection
      if(td.gt.tf)td=tf         !if cluster forms within td

      call read()
       
!     calculate the total number of mergers at a (look-back) time <t
      open(unit=77,file="Nmergers.txt",action="WRITE")
      write(77,*)"time[Gyr], N_in(<t), N_out(<t)"
      
      eccen=1.d5
      dt=0.05                   !output stepsize
      t=0.d0      
      i=0
      do 
         i=i+1                 
         td0=tf-t
         t=t+dt*1.d9
         if(td0.lt.0)exit 
         
!     in-cluster mergers
         X0=0.
         X1=td0                
         call dqag (fun,X0,X1, epsabs, epsrel, key, Fin,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )
         
!     ejected binaries
         X0=0.
         X1=td0                
         call dqag (funej,X0,X1, epsabs, epsrel, key, Fej1,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )
         
         X0=td0
         X1=tf        
         call dqag (funej2,X0,X1,epsabs, epsrel, key, Fej2,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )
         
         write(77,"(F14.2,2F14.7)")td0/1.d9,Fin,Fej1+Fej2
         end do
      
      
!     Cumulative eccentricity distribution of binaries that merge at t<td
      open(unit=99,file='Necc.txt',action="WRITE")
      write(99,*)"e, N(<e)"     

      X0=0.d0
      X1=td 
      call dqag (fun,X0,X1,epsabs, epsrel, key, Fin,abserr, 
     &     neval, ier, limit, lenw, last, iwork, work )
      X0=0.d0
      X1=td                
      call dqag (funej,X0,X1, epsabs, epsrel, key, Fej1,abserr, 
     &     neval, ier, limit, lenw, last, iwork, work )
      X0=td
      X1=tf
      call dqag (funej2,X0,X1,epsabs, epsrel, key, Fej2,abserr, 
     &     neval, ier, limit, lenw, last, iwork, work )            
      Ntot=Fin+Fej1+Fej2        !total number of mergers within td
      
      step=0.01
      dec=0.
      do while(Ntot.gt.0)
         dec=dec+step
         eccen=10**(-10.+dec)

!     in-cluster mergers
         X0=0.
         X1=td                
         call dqag (fun,X0,X1,epsabs, epsrel, key, Fin,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )

!     ejected binaries
         X0=0.d0
         X1=td               
         call dqag (funej,X0,X1, epsabs, epsrel, key, Fej1,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )        
         X0=td
         X1=tf
         call dqag (funej2,X0,X1,epsabs, epsrel, key, Fej2,abserr, 
     &        neval, ier, limit, lenw, last, iwork, work )               
         
         if(eccen.ge.1.d0)exit
         write(99,"(ES14.7,F14.7)")eccen,(Fej1+Fej2+Fin)/Ntot
      enddo
      
      
      end program  BHBdynamics

            
      
