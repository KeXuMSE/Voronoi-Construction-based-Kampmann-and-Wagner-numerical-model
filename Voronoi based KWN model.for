      program kwn !Ni-7.5at%Al-8.5at%Cr at 600C      
C...define parameters    
      implicit double precision (a-h, o-z),integer(i-n)
C...parameters about growth, rate nucleation rate etc.      
      double precision Gmgp,Q,DGM,DGV,GNV,DM,RC,REFF
      double precision nig,alg,crg,nigp,algp,crgp,nige,alge,crge
      double precision gnimu,galmu,gcrmu
      double precision rate,zeld,belta
C...number of new gamma prime and total gamma prime in 100^3 nm cell 
      double precision ngpcell,tgpcell,tgpsum
C...iteration times
      integer id,jd,kd,ld,md
      integer ka,la
      integer idi
C...sum and average of effective diffusion distance
      double precision EDSUM,EDAVG
C...result of evolution N f NiRi R alau algp crau crgp     
      double precision RES(0:1000000,1:8)   !!!change time interval
C...EVR array i.e. radius array in 100nm^3 cell
      double precision EVR(1,100000)        !!!change array size
C...EVO array i.e. X Y Z R sigmar(effective diffusion distance)
      double precision EVO(1:5,100000)      !!!change array size
C...define the X,Y,Z,DIST,COMB,COMS to calculate the xyz array
      double precision X,Y,Z,COMB,COMS
      double precision DIST(1,100000)       !!!change array size
C...define the Input and Output array of matlab
      double precision ainp(100000,3),bout(100000,1)    !!!notice column and row   
C...critical nucleus compositions
      double precision cnal,cncr
C...define parameters about gfile      
      parameter(nwg=800000,nwp=5000000)
      dimension iwsg(nwg),iwse(nwp)
      character*256 tcpath,tmppath
      integer iwsg,iwse
      character*32 gfile
      logical sg2err 
C...define parameters about matlab interface      
      integer*8,external::engOpen,engClose,mxCreateDoubleMatrix
      integer*8,external::mxGetPr,engPutVariable,engGetVariable
      integer*8,external::engEvalString
      integer*8 ep,T
      integer status           

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
      tgpcell=0D0   !!!change initial number of gamma prime in 100^3 nm cell 
C...al and cr in full austenite at time 0  !!!change intial label and value    
      RES(0,5)=0.075D0
      RES(0,7)=0.085D0
C...al and cr in gamma prime at time 0
      RES(0,6)=0D0
      RES(0,8)=0D0
C...composition in critical nucleus
      cnal=0.168D0
      cncr=0.081D0      
C...N,f,sum(Ri),Ravg at time 0
      RES(0,1)=0D0   !!!change initial number density
      RES(0,2)=0D0
      RES(0,3)=0D0
      RES(0,4)=0D0   !!!change initial average size
C...define the initial value X Y Z R for calculating the DIST
      EVO(1,0)=0.0D0  !X
      EVO(2,0)=0.0D0  !Y
      EVO(3,0)=0.0D0  !Z
      EVO(4,0)=0.0D0  !R
      EVO(5,0)=0.0D0  !sigma r effective diffusion distance



C...comments
      write(*,1)
1     format('Ni-Al-Cr at 600C evolution')      
C...open print file
      open(2,file='C:\Users\Public\Documents\Thermo-Calc\2021a\SDK\TQ\
     &Windows\fortran\KWN\kwnresult.txt')
      open(3,file='C:\Users\Public\Documents\Thermo-Calc\2021a\SDK\TQ\
     &Windows\fortran\KWN\kwnnr.txt')
      
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
C...open matlab interface and set working path      
C...open matlab
      ep=engOpen('')
      if (ep == 0) then
          write(*,*) 'can''t start matlab engine'
          stop
      end if
C...set path
      if (engEvalString(ep," addpath('C:\Program Files\MATLAB\R2017b\
     &toolbox\mpt\mpt')") /= 0) then
          write(*,*) 'engEvalString failed'
          stop
      end if   
      
C...evolution(time interval changes)      
      do 10 i=1,1000000                        !!!change time interval 1000000s
          if (i.LE.1000000D0) then             !!!change time interval 1000000s
              tevo=1.0D0*i                     !!!change time interval 1s     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C...get the EVR array in 100nm^3 cell
C...get the nulceation rate
              ngpcell=ANINT(1.0D-21*DM)   !!! new added number   !!!change time interval 1s
              tgpcell=tgpcell+ngpcell
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C...get the EVR array             
              do 11 j=1,tgpcell
                  if (j .LE. tgpcell-ngpcell) then
                      EVR(1,j)=EVR(1,j)
                  else
                      EVR(1,j)=REFF
                  end if
11            continue                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...get the random coordination XYZ and R array (i.e. EVO array)
              do 12  ka=1,tgpcell
C...get the radius
                  EVO(4,ka)=EVR(1,ka)
C...generate random coordiantion
121               call random_number(X)
                  call random_number(Y)
                  call random_number(Z)
C...check whether every former point satisfies D/2 .GL. Rbigger & D/2 .LT. Rsmaller
                  do 122 la=1,ka
                      DIST(1,la)=sqrt((1D-7*X-EVO(1,la-1))**2
     &+(1D-7*Y-EVO(2,la-1))**2+(1D-7*Z-EVO(3,la-1))**2)   
                      COMB=0.5D0*DIST(1,la)-EVO(4,la)
                      COMS=0.5D0*DIST(1,la)-EVO(4,la-1)
                      if(COMB .GT. 0.0D0 .AND. COMS .GT. 0.0D0) then
                          go to 122
                      else 
                          go to 121                  
                      end if
122               continue  
C...make sure the coordination and get the array ainp
                  EVO(1,ka)=1D-7*X
                  EVO(2,ka)=1D-7*Y
                  EVO(3,ka)=1D-7*Z  
                  ainp(ka,1)=1D-7*X
                  ainp(ka,2)=1D-7*Y
                  ainp(ka,3)=1D-7*Z
12            continue                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC              
C...set the rest of ainp as 0
              do 13 l=tgpcell+1,100000   !!!change size of ainp array
                  ainp(l,1)=0.0D0
                  ainp(l,2)=0.0D0
                  ainp(l,3)=0.0D0
13            continue              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...calculate effective diffusion distance by voronoi in matlab             
C...get the ainp array (x,y,z) in matlab from Fortran
              T=mxCreateDoubleMatrix(100000,3,0)     !!!change size of ainp array
              call mxCopyReal8ToPtr(ainp,mxGetPr(T),100000*3)   !!!change size of ainp array
              status=engPutVariable(ep,'T',T)
              if (status /= 0) then
                  write(*,*) 'engPutVariable failed'
                  stop
              end if
C...give the element, which is not equeal to 0 in array A, to array B
              if (engEvalString(ep,"A=1E9*T;for m=1:100000;if(A(m,1)>0);        !!!change size of ainp array
     &B(m,1)=A(m,1);end;if(A(m,2)>0);B(m,2)=A(m,2);end;
     &if(A(m,3)>0);B(m,3)=A(m,3);end;end") /= 0) then
                  write(*,*) 'engEvalString failed'
                  stop
              end if 
C...draw the voronoi construction by B array 
              if (engEvalString(ep,"Options.plot=1;v=[0 0 0;0 0 100;
     &0 100 0;0 100 100;100 0 0;100 0 100;100 100 0;100 100 100];
     &P=polytope(v);Options.pbound=P;
     &Pn=mpt_voronoi(B,Options)") /= 0) then
                  write(*,*) 'engEvalString failed'
                  stop
              end if   
C...calculate effective diffusion distance     
              if (engEvalString(ep,"for n=1:size(B,1);V=extreme(Pn(n));
     &for o=1:size(V,1);
     &D(1,o)=norm([B(n,1),B(n,2),B(n,3)]-[V(o,1),V(o,2),V(o,3)]);
     &end;E(n,1)=sum(D)/size(V,1);clear D;clear V;end") /= 0) then
                  write(*,*) 'engEvalString failed'
                  stop
              end if  
C...calculate E and bout array    
              if (engEvalString(ep,"for p=1:100000;if(p>size(B,1));    !!!change size of ainp array
     &E(p,1)=0;end;end") /= 0) then
                  write(*,*) 'engEvalString failed'
                  stop
              end if                      
C...get the result in Fortran from matlab
              status = engGetVariable(ep,"E")  
              call mxCopyPtrToReal8(mxGetPr(status),bout,100000*1)     !!!change size of ainp array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
C...get the EVO(5,id) array effective diffusion distance            
              do 14 id=1,tgpcell
                  EVO(5,id)=1.0D-9*bout(id,1)-EVO(4,id) 
14            continue
C...calculate average of diffusion distance
              EDSUM=0.0D0
              do 141 idi=1,tgpcell
                  EDSUM=EDSUM+EVO(5,idi)    
141           continue
              EDAVG=EDSUM/tgpcell
C...clear all
              if (engEvalString(ep,"clc;clear all;close all") /= 0) then
                  write(*,*) 'engEvalString failed'
                  stop
              end if 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C...update the EVO by growth rate equation (EVO(4,jd))      
              do 15 jd=1,tgpcell
                  EVO(4,jd)=EVO(4,jd)+
     &1.0D0*(DGM-2*SFE*VM/EVO(4,jd))/rate/EVO(5,jd) 
C...if radius is less than 0, set EVO(4,jd) as 0                  
                  if (EVO(4,jd).LT.1.0D-10) then
                      EVO(4,jd)=0D0
                  end if   
15            continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C... define initial values of RES
              RES(i,2)=0
              RES(i,3)=0      
C...calculate the RES array      
              do 16 kd=1,tgpcell
C...sum f (not unity)
                  RES(i,2)=RES(i,2)+4*PI*EVO(4,kd)*EVO(4,kd)*EVO(4,kd)/3
C...sum (Ri)
                  RES(i,3)=RES(i,3)+EVO(4,kd)
16            continue
C...sum N in 1m^3
              RES(i,1)=1.0D21*tgpcell     
C...Ravg     
              RES(i,4)=RES(i,3)/tgpcell
C...sum f in 1m3
              RES(i,2)=1.0D21*RES(i,2)
C...aveal in gamma prime
              RES(i,6)=RES(i-1,6)*RES(i-1,2)/RES(i,2)+
     &4*PI*cnal*(REFF**3)*1.0D0*DM/RES(i,2)/3+
     &algp*(RES(i,2)-RES(i-1,2)-4*PI*(REFF**3)*1.0D0*DM/3)/RES(i,2)
C...aual
              RES(i,5)=(RES(0,5)-RES(i,6)*RES(i,2))/(1-RES(i,2))
C...avecr in gamma prime
              RES(i,8)=RES(i-1,8)*RES(i-1,2)/RES(i,2)+
     &4*PI*cncr*(REFF**3)*1.0D0*DM/RES(i,2)/3+         
     &crgp*(RES(i,2)-RES(i-1,2)-4*PI*(REFF**3)*1.0D0*DM/3)/RES(i,2)
C...aucr
              RES(i,7)=(RES(0,7)-RES(i,8)*RES(i,2))/(1-RES(i,2))
              print 22,tevo,RES(i,1),RES(i,2),RES(i,4),REFF,DM,DGM,GNV,
     &RES(i,5),RES(i,7),RES(i,6),RES(i,8),EDAVG  
              write(2,22),tevo,RES(i,1),RES(i,2),RES(i,4),REFF,DM,DGM,
     &GNV,RES(i,5),RES(i,7),RES(i,6),RES(i,8),EDAVG     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...update the EVR array and tgpcell
C...define initial value
              tgpsum=0.0D0
              do 17 ld=1,tgpcell
                  if (EVO(4,ld) .NE. 0D0) then
                      tgpsum=tgpsum+1
                      EVR(1,tgpsum)=EVO(4,ld)                      
                  end if                  
17            continue
C...get the new number in 100nm^3
              tgpcell=tgpsum 
C...print EVR
              write(3,33),tevo,(EVR(1,md),md=1,tgpcell)
33            format(100001ES11.3)              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...end if of iteration 10
          end if            
10    continue 
            
            
C...print content     
22    format(13ES11.3) 
                 
      pause
      end program kwn
      
      


      




