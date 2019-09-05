!!!!! evolvegcbh: simple model for cluster evolution
      subroutine simple(Nmax,rh_i,Mcl_i,fbh0)
      USE modelPar
      USE inputPar
      implicit double precision (a-z)
      external evolve
      integer IR,nq,l,i,j,merge,im,im2,k,h,nbh,kick,Nmax
      parameter (nq=3)
      dimension y(nq),tol(nq)
      dimension radius(Nmax),Mass(Nmax),mass_max(Nmax)
      dimension time(Nmax),massbh(Nmax)
      data tol/nq*1d-18/
      common/cnst/pi 
      common/tbe/tb,psi
      psi=5.d0                  !relaxation parameter
      open (unit=88,file="cluster.txt",action="WRITE") !file with: time, M_bh, M_cl, r_h, m_max
     
      
      DT0=1.                 !initial timestep
 500  y(1)=rh_i
      y(2)=Mcl_i
      y(3)=Mcl_i*fbh0      
      call TBalanced(y(3),y(2),y(1),tb)      
      
      time=0.d0
      massbh=0.d0
      mass=0.d0
      radius=0.d0
      radius(1)=rh_i
      massbh(1)=Mcl_i*fbh0
      mass(1)=Mcl_i
      mass_max(1)=funm(1.d0)
      
      DT=tb
      t=tb                    
      tmax=tf/1.e9  
      i=1      
      do   while(t.lt.tmax)
         call RK78(IR,t,dt,y,3,tol,evolve)
         if(y(1).ne.y(1))then
            dt0=dt0/10.
            goto 500
         end if         
         if(y(3).lt.0)y(3)=0.d0
         i=i+1
         time(i)=t         
         massbh(i)=y(3)
         mass(i)=y(2)
         radius(i)=y(1)
         P=massbh(i)/(Mcl_i*fbh0)
         mass_max(i)=funm(P)    !deterime the mass of the BHs
      end do

      do j=1,i
         write(88,*)time(j),massbh(j),mass(j),radius(j),mass_max(j)
      end do
     
      
      close(88)
      end


      subroutine evolve(t,y,yp)
      USE modelPar
      USE inputPar
      implicit double precision (a-z)
      integer n,k,h
      parameter  (n=3)
      dimension y(n),yp(n)
      common/cnst/pi
      common/tbe/tb,psi      
    
!     only relaxation dirven evolution included
      trh=2.06d5*sqrt(y(2)*y(1)**3)/psi/m_mean/1.d9 !relaxation time            

!     expansion and BHmass loss due to relaxation      
      yp(3)=-beta*y(2)/trh      !mass loss of BHs
      yp(2)=yp(3)               !change of total mass 
      yp(1)=(csi/trh+2.d0*yp(2)/y(2))*y(1) !radius expansion
      if(t.lt.tb)yp=0.
      end
      
