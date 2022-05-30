      program kwn !Ni-7.5at%Al-8.5at%Cr at 600C      
C...define parameters    
      implicit double precision (a-h, o-z),integer(i-n)
C...parameters about growth rate & nucleation rate etc.      
      double precision Gmgp,Q,DGM,DGV,GNV,DM,RC,REFF
      double precision nig,alg,crg,nigp,algp,crgp,nige,alge,crge
      double precision gnimu,galmu,gcrmu
      double precision rate,zeld,belta
C...result of evolution N f NiRi R alau algp crau crgp     
      double precision RES(0:3686400,1:8)   !!!change time interval
C...evolution array N:EVO(i,1) R:EVO(i,2)   
      double precision EVO(1:3686400,1:2)   !!!change time interval
C...critical nucleus compositions
      double precision cnal,cncr
C...define parameters about gfile      
      parameter(nwg=800000,nwp=5000000)
      dimension iwsg(nwg),iwse(nwp)
      character*256 tcpath,tmppath
      integer iwsg,iwse
      character*32 gfile
      logical sg2err       
C...define constants
      QNI=288310D0                    !J/mol
      QAL=281590D0                    !J/mol
      QCR=281150D0                    !J/mol
      DNI=9.64D-22                    !M**2/s Ni in Ni at 873K
      DAl=8.31D-21                    !M**2/s Al in Ni at 873K
      DCR=4.39D-21                    !M**2/s Cr in Ni at 873K      
      PI=3.141592654D0
      SFE=0.0285D0                    !J/M**2   surface energy  !!!change
      BOL=1.3806505D-23				!JOL/K     
      DMS=3.0D27                      !density of nucleation site  !!!change
 	VM=6.59D-6			            !M**3/mol   molar volume of gamma prime     
      Va=1.09D-29                     !M**3   volume per atom in the matrix
      alattice=3.4D-10                !M lattice constant of gamma prime       
      Temp=873.15D0                   !K
C...define initial values
      tevo=0D0
C...al and cr in full austenite at time 0  !!!change intial value   
      RES(0,5)=0.075D0
      RES(0,7)=0.085D0
C...al and cr in gamma prime at time 0
      RES(0,6)=0D0
      RES(0,8)=0D0
C...N,f,NiRi,R at time 0
      RES(0,1)=0D0
      RES(0,2)=0D0
      RES(0,3)=0D0
      RES(0,4)=0D0            
C...composition in critical nucleus
      cnal=0.168D0
      cncr=0.081D0      
C...comments
      write(*,1)
1     format('Ni-Al-Cr at 600C evolution')      
C...open print file
      open(2,file='C:\Users\Public\Documents\Thermo-Calc\2021a\SDK\TQ\
     &Windows\fortran\KWN\kwnresult.txt')
C...open gfile   
C...gfile
      gfile='NIALCR.GES5'                 
C...initiate the workspace
      tcpath=' '
      tmppath=' '
      call tqini3(tcpath,tmppath,nwg,nwp,iwsg,iwse)
C...read the thermodynamic data file which was created by using
C...the GES module inside the Thermo-Calc software package
      call tqrfil(gfile,iwsg,iwse)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C...evolution      
      do 10 i=1,3686400                        !!!change time interval 3686400s
          if (i.LE.3686400D0) then             !!!change time interval 3686400s
              tevo=1.0D0*i                     !!!change time interval 1s     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...get new nucleus density reff and etc.(subrountine)              
C...get component indexes      
              call tqgsci(indexni,'NI',iwsg,iwse)
              call tqgsci(indexal,'AL',iwsg,iwse)
              call tqgsci(indexcr,'CR',iwsg,iwse)
C...get phase indexes      
              call tqgpi(indexg,'FCC_L12',iwsg,iwse)
              call tqgpi(indexgp,'FCC_L12#2',iwsg,iwse)
C...set T/K
              call tqssu('T','K',iwsg,iwse)      
C...set condition of equilibrium 
              call tqsetc('T',-1,-1,Temp,icont,iwsg,iwse)
              call tqsetc('N',-1,-1,1.00D0,iconn,iwsg,iwse)
              call tqsetc('P',-1,-1,101325.0D0,iconp,iwsg,iwse)
              call tqsetc('X',-1,indexal,RES(i-1,5),iconal,iwsg,iwse)   !!!change initial value
              call tqsetc('X',-1,indexcr,RES(i-1,7),iconcr,iwsg,iwse)   !!!change initial value
C...calculate equilibrium
              call tqce(' ',0,0,0.0D0, iwsg,iwse)
              if(sg2err(ierr)) print *, 'Calculation failed!'                
C...get mole fraction of gamma prime      
              call tqget1('X',indexgp,indexni,nigp,iwsg,iwse)
              call tqget1('X',indexgp,indexal,algp,iwsg,iwse)
              call tqget1('X',indexgp,indexcr,crgp,iwsg,iwse)
C...get mole fraction of gamma after precipitation
              call tqget1('X',indexg,indexni,nige,iwsg,iwse)
              call tqget1('X',indexg,indexal,alge,iwsg,iwse)
              call tqget1('X',indexg,indexcr,crge,iwsg,iwse)               
C...get Gm of fcc_l12
              Gmgp=tqggm(indexgp,iwsg,iwse)
C...calculate deltaG* and Qeff
              call tqcsp(indexg,'ENTERED',1.0D0,iwsg,iwse)
              call tqcsp(indexgp,'SUSPENDED',0.0D0,iwsg,iwse)
C...calculate equilibrium
              call tqce(' ',0,0,0.0D0,iwsg,iwse)
              if(sg2err(ierr)) print *, 'Calculation failed!'       
C...get mole fraction of elements in fcc_a1
              call tqget1('X',-1,indexni,nig,iwsg,iwse)
              call tqget1('X',-1,indexal,alg,iwsg,iwse)
              call tqget1('X',-1,indexcr,crg,iwsg,iwse)      
C...get Qeff
              Q=QNI+alg*QAL+crg*QCR      
C...get chemical potential of elements in fcc_a1
              call tqget1('MU',-1,indexni,gnimu,iwsg,iwse)
              call tqget1('MU',-1,indexal,galmu,iwsg,iwse)
              call tqget1('MU',-1,indexcr,gcrmu,iwsg,iwse)     
C...calculate deltaGm deltaGv
              DGM=nigp*gnimu+algp*galmu+crgp*gcrmu-Gmgp
              DGV=DGM/VM      
C...get deltaG*     
              GNV=16*PI*(SFE**3)/DGV/DGV/3
C...get critical nucleation radium and reff
              RC=2*SFE/DGV
              REFF=RC+0.5*((BOL*Temp/PI/SFE)**0.5)           
C...get nucleation rate      
              zeld=Va*DGV*DGV/8/PI/((SFE*SFE*SFE*BOL*Temp)**0.5)     !zeldovich factor
C...rate of attachment of solute atoms    
              belta=4*PI*RC*RC/alattice**4/
     &(1/DNI/nig+1/DAL/alg+1/DCR/crg)  
C...nucleation rate
              DM=DMS*(1-RES(i-1,2))*zeld*belta*exp(-GNV/BOL/Temp)      
C...rate of growth    
              rate=8.314D0*Temp*((nigp-nige)**2/nige/DNI+(algp-alge)**2/
     &alge/DAL+(crgp-crge)**2/crge/DCR)             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...get new nucleus density and reff              
              EVO(i,1)=1.0D0*DM                     !!!change time interval
              EVO(i,2)=REFF              
C...calculate the EVO array              
              do 101 j=1,i
                  if (EVO(j,1).NE.0D0) then
                      EVO(j,2)=EVO(j,2)+
     &1.0D0*(DGM-2*SFE*VM/EVO(j,2))/rate/EVO(j,2)   !!!change time interval
C...if particle size is less than 1D-10, or particle is less than 0D0, then erase this particle group                   
                      if ((EVO(j,1).LT.0.0D0).OR.
     &(EVO(j,2).LT.1.0D-10)) then
                          EVO(j,1)=0D0
                          EVO(j,2)=0D0
                      end if
                  end if
101           continue              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C... define initial values of RES
              RES(i,1)=0
              RES(i,2)=0
              RES(i,3)=0                
C...calculate the RES array      
              do 102 k=1,i
C...sumN                
                  RES(i,1)=RES(i,1)+EVO(k,1)
C...sumf
                  RES(i,2)=RES(i,2)+4*PI*EVO(k,2)*EVO(k,2)*EVO(k,2)*
     &EVO(k,1)/3
C...NiRi
                  RES(i,3)=RES(i,3)+EVO(k,1)*EVO(k,2)
102           continue
C...aveR              
              RES(i,4)=RES(i,3)/RES(i,1)
C...aveal in gamma prime            
              RES(i,6)=RES(i-1,6)*RES(i-1,2)/RES(i,2)+
     &4*PI*cnal*(REFF**3)*1.0D0*DM/RES(i,2)/3+                          !!!change time interval
     &algp*(RES(i,2)-RES(i-1,2)-4*PI*(REFF**3)*1.0D0*DM/3)/RES(i,2)     !!!change time interval         
C...aual
              RES(i,5)=(RES(0,5)-RES(i,6)*RES(i,2))/(1-RES(i,2))
C...avecr in gamma prime
              RES(i,8)=RES(i-1,8)*RES(i-1,2)/RES(i,2)+
     &4*PI*cncr*(REFF**3)*1.0D0*DM/RES(i,2)/3+                          !!!change time interval
     &crgp*(RES(i,2)-RES(i-1,2)-4*PI*(REFF**3)*1.0D0*DM/3)/RES(i,2)     !!!change time interval
C...aucr
              RES(i,7)=(RES(0,7)-RES(i,8)*RES(i,2))/(1-RES(i,2)) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                            
C...print result              
              print 22,tevo,RES(i,1),RES(i,2),RES(i,4),REFF,DM,DGM,GNV,
     &RES(i,5),RES(i,7),RES(i,6),RES(i,8)  
              write(2,22),tevo,RES(i,1),RES(i,2),RES(i,4),REFF,DM,DGM,
     &GNV,RES(i,5),RES(i,7),RES(i,6),RES(i,8)     
C...end if of iteration 10
          end if            
10    continue             
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
C...print content     
22    format(12ES11.3) 
                 
      pause
      end program kwn
      
      


      




