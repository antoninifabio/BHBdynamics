!     COMMON block module2 
      MODULE inputPar
      integer input,IMF
      real*8 freq,tf,a
      SAVE input,IMF,freq,tf,a
      END MODULE inputPar
       
      MODULE commonV      
      integer, parameter:: Nmax=100000
      integer cont,cont2
      real*8,dimension(Nmax):: tabt,tabr,tabm,tabmbh
      real*8,dimension(Nmax):: tabmbh_p,tabmax
      real*8,dimension(Nmax):: c1,c2,c3
      real*8,dimension(Nmax):: c4,c5,c6
      real*8,dimension(Nmax):: c7,c8,c9
      real*8,dimension(Nmax):: c10,c11,c12
      real,dimension(Nmax):: w
      integer,dimension(Nmax):: INDEX
      SAVE tabt,tabm,tabr,tabmbh,tabmbh_p,tabmax,cont,cont2
      SAVE c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,INDEX,w
      END MODULE commonV
      
      MODULE modelPar      
      real*8 csi,m_mean,tcc,a1,beta,dE     
      SAVE csi,m_mean,tcc,a1,beta,dE     
      END MODULE modelPar
