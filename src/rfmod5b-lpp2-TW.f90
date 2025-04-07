

MODULE numberformat
!--------------------
!-------
   IMPLICIT none
   INTEGER, PARAMETER::rk=kind(1d0)
END MODULE numberformat
!--------------------


MODULE constants
!--------------------
!-------
! Contains constants for calculations
!-------
   USE numberformat
   IMPLICIT none
   REAL(rk),PARAMETER::c=2.9979e10_rk ! is the speed of light [cm/s]
   REAL(rk),PARAMETER::pi=3.141592654_rk
   REAL(rk),PARAMETER::mu0=4.0_rk*pi*1e-9_rk !is the magnetic conductivity
   REAL(rk),PARAMETER::eps0=8.85e-14_rk !is the permittivity of free space
   REAL(rk),PARAMETER::ec=4.8088e-10_rk !is the electron charge
   REAL(rk),PARAMETER::em=9.1094e-28_rk !is the electron mass [g]
   REAL(rk),PARAMETER::pm=1.6726e-24_rk !is the ion mass [g]
!REAL(rk),PARAMETER::pm=2.0_rk*1.6726e-24_rk !is the double ion mass [g]
   REAL(rk),PARAMETER::CL=13 !Coulomb logarithm
   REAL(rk),PARAMETER::CONBOL=1.602e-12_rk !is the Boltzman constant [erg/eV]
!REAL(rk),PARAMETER::SVEN=1e-7_rk !is the rate of the electron-neutral collisions
!REAL(rk),PARAMETER::SVIN=1e-8_rk !is the rate of the ion-neutral collisions

END MODULE constants
!--------------------



!--------------------
MODULE eqpar
!--------------------
!-------
!numeq - number of equations,
!numfem - number of finite elements per node,
!nfm3,nfm3 - number of finite elements at the segment
!-------
   IMPLICIT none
   INTEGER, PARAMETER::	numeq=3,  &
      numfem=2, &
      numline=numeq*numfem, &
      nfm3=2*numfem
END MODULE eqpar
!--------------------


!--------------------
MODULE gaussdata
!--------------------
!-------
!Contains constants for Gauss quadrature formula
!-------
   USE numberformat
   IMPLICIT none
   INTEGER, PARAMETER::kgmax=7
   REAL(rk),PARAMETER,PRIVATE::w1=0.417959183673469_rk, &
      w2=0.381830050505119_rk, &
      w3=0.279705391489277_rk, &
      w4=0.129484966168870_rk, &
      x1=0.405845151377397_rk, &
      x2=0.741531185599394_rk, &
      x3=0.949107912342559_rk

   REAL(rk), DIMENSION(kgmax),PARAMETER::weigh=(/w4/2.0_rk, &
      w3/2.0_rk, &
      w2/2.0_rk, &
      w1/2.0_rk, &
      w2/2.0_rk, &
      w3/2.0_rk, &
      w4/2.0_rk/), &
      ksi=(/(-x3+1.0_rk)/2.0_rk, &
      (-x2+1.0_rk)/2.0_rk, &
      (-x1+1.0_rk)/2.0_rk, &
      .5_rk,          &
      (x1+1.0_rk)/2.0_rk,  &
      (x2+1.0_rk)/2.0_rk,  &
      (x3+1.0_rk)/2.0_rk/)
END MODULE gaussdata
!--------------------


!--------------------
MODULE rfdata
!--------------------
!-------
!Contains calculation parameters
!-------
   USE numberformat
   IMPLICIT none
   INTEGER::numptinit,nant,mmin,mmax,nmin,nmax !initial number of mesh points, number of antennas,
   !minimum and maximum azimuthal and axial mode number
   INTEGER::numpt                           !actual number of mesh points
   REAL(rk)::a !is the radius of metallic wall [cm]
   REAL(rk)::rtor !is the major radius of equivalent torus [cm]
   REAL(rk)::rstar !is the equivalent radius in vertical direction [cm]
   REAL(rk)::w0 ! is the RF frequency [s-1]
   REAL(rk)::radac,width,weight !parameters for mesh accumulation(mesh accumulation radius [cm],
   !width of the region [cm], weight (0...1)).

   REAL(rk)::khmax, khmin,epspmin !K_perp*h product (how rare is mesh), theshold for patching Epsylon_perp

END MODULE rfdata
!--------------------

MODULE antdata
!--------------------
!-------
!Contains calculation parameters
!-------
   USE numberformat
   IMPLICIT none
   INTEGER, PARAMETER::maxnumant=20
   INTEGER::nanten,nanthfs ! number of antennas, number of hfs antennas
   INTEGER,DIMENSION(maxnumant)::antnode !mesh node number corresponding to each antenna
   REAL(rk),DIMENSION(maxnumant)::radant !antenna radii (the same as rad(antnode())) [cm]
   REAL(rk),DIMENSION(maxnumant)::thetaant !theta position of antenna
   REAL(rk),DIMENSION(maxnumant)::phiant !azimutal position of antenna
   REAL(rk),DIMENSION(maxnumant)::delphiant !azimutal size of antenna stap
   REAL(rk),DIMENSION(maxnumant)::delthetaant !theta size of antenna
   REAL(rk)::sldown !slowing down coefficient of the TEM wave along the strap
   COMPLEX(rk),DIMENSION(maxnumant)::phase ! antenna strap current multiplier
   INTEGER,DIMENSION(maxnumant)::ordering !antenna numbers in order at which radant increase

END MODULE antdata



!--------------------
SUBROUTINE RFpower (hfrq,zsort, msort, nsort, rad1, numpt1, density, magnf, cfreq, powin,alr, powr,lastcall)
!--------------------
!Zsort(nsort) �  ion charges in units of proton charge;
!Msort(nsort) � ion mass number;
!Nsort -             number of sorts of ions;
!Rad1(npoints)-mesh radial points (-a:a) [cm];
!Numpt1 -         number of mesh points;
!density(npoints,nsort)- ion density distributions [cm-3];
!Magnf(npoints) � magnetic field radial profile [G];
!cfreq(npoints,nsort+1)- collision frequencies [c-1];
!powin �       input power [W];
!alr �       antenna loading resistance [Ohm] (output);
!powr(npoints,nsort+2)- RF power profile delivered to each ion sort, electrons, and full RF power including power due to penaltizing [W/cm3](output);
!Lastcall � .false. if there will be more calls of the subroutine and .true. if the call is last.
!-------
   USE numberformat
   USE eqpar
   USE rfdata
   USE constants
   IMPLICIT none
   INTEGER::i,k
   REAL(rk), DIMENSION(:),ALLOCATABLE,SAVE::rad
   REAL(rk), DIMENSION(:),ALLOCATABLE::rad7
   REAL(rk),INTENT(IN)::hfrq
   INTEGER,INTENT(IN)::numpt1,nsort
   COMPLEX(rk),DIMENSION(:,:),ALLOCATABLE,SAVE::sysmatr,amatr1p,amatr1m,amatr0,bmatr
   COMPLEX(rk),DIMENSION(:),ALLOCATABLE,SAVE::yy
   COMPLEX(rk),DIMENSION(:,:),ALLOCATABLE,SAVE::epsparr,epsgarr,epslarr
   COMPLEX(rk),DIMENSION(:,:),ALLOCATABLE::epsparr7,epslarr7
   COMPLEX(rk),DIMENSION(:,:,:),ALLOCATABLE,SAVE::rhsarr
   COMPLEX(rk)::epsp,epsg,epsl
   REAL(rk),ALLOCATABLE,DIMENSION(:,:),SAVE::p
   REAL(rk),ALLOCATABLE,DIMENSION(:,:),SAVE::ptmp
   REAL(rk),ALLOCATABLE,DIMENSION(:,:,:),SAVE::powm
   REAL(rk)::tp
   REAL(rk),DIMENSION(numpt1)::rad1
   REAL(rk),DIMENSION(2*numpt1,nsort)::density
   REAL(rk),DIMENSION(numpt1,nsort+1)::cfreq
   REAL(rk),DIMENSION(numpt1)::magnf
   REAL(rk),DIMENSION(numpt1,nsort+2)::powr
   REAL(rk),INTENT(IN)::powin
   REAL(rk),INTENT(OUT)::alr
   LOGICAL::lastcall
   LOGICAL,SAVE::firstcall=.TRUE.
   INTEGER::m,n
   INTEGER,DIMENSION(nsort)::zsort,msort
   INTEGER::mmaxn
   REAL(rk),ALLOCATABLE,DIMENSION(:),SAVE::kper2
   REAL(rk)::kpermax,coeffi
   INTEGER::fidwr																						! TW20190225

   interface																									! TW20190225
      integer function newunit(unit)													! TW20190225
         integer, intent(out), optional :: unit								! TW20190225
      end function newunit														  			! TW20190225
   end interface																							! TW20190225

   IF(firstcall) THEN
!
! Reading data from files
!
      CALL inrf
      CALL inant
      firstcall=.FALSE.
      w0=hfrq*1e6_rk*2.0_rk*pi
      numpt=numptinit
      a=rad1(numpt1)

!-------
!Memory allocation
!-------

      ALLOCATE(rad(numpt))
      ALLOCATE(amatr1p(3*numline,numline*numpt))
      ALLOCATE(amatr1m(3*numline,numline*numpt))
      ALLOCATE(amatr0(3*numline,numline*numpt))
      ALLOCATE(bmatr(3*numline,numline*numpt))
      ALLOCATE(sysmatr(3*numline,numline*numpt))
      ALLOCATE(yy(numline*numpt))
      ALLOCATE(epsparr(numpt,nsort+2),epsgarr(numpt,nsort+2),epslarr(numpt,nsort+2))
      ALLOCATE(rhsarr(numline*numpt,mmin:mmax,nmin:nmax))
      ALLOCATE(p(numpt,nsort+2))
      ALLOCATE(ptmp(numpt,nsort+2))
      ALLOCATE(powm(numpt,nsort+2,mmin:mmax))
      ALLOCATE(kper2(numpt))

!-------
!Creating the mesh
!-------
      CALL mesh(rad)
!-------
!Preparing the right-hand-side
!-------

      DO n=nmin,nmax
         DO m=mmin,mmax

            CALL rhsmultistrap(rad,yy,m,n)

            rhsarr(:,m,n)=yy

         END DO
      END DO


   END IF

!-------
!Power density array initiation
!-------
   p=0.0_rk

!-------
!Epsylon-dependent part of the matrix
!-------

   CALL projeqeps(rad,rad1,numpt1,bmatr,density,cfreq,magnf,zsort,msort,nsort)
!-------
!Matrix for power calculations
!-------
   DO i=1,numpt
      DO k=1,nsort

         CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,k,k,.FALSE.)

         epsparr(i,k)=epsp
         epsgarr(i,k)=epsg
         epslarr(i,k)=epsl

      END DO

      CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,1,0,.TRUE.)

      epsparr(i,nsort+1)=epsp
      epsgarr(i,nsort+1)=epsg
      epslarr(i,nsort+1)=epsl

      CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,1,nsort,.TRUE.)

      epsparr(i,nsort+2)=epsp
      epsgarr(i,nsort+2)=epsg
      epslarr(i,nsort+2)=epsl

   END DO
!
! Block for determinig number of mesh points
!

   DO i=1,numpt
      kper2(i)=ABS(-epslarr(i,nsort+2)/epsparr(i,nsort+2)*(max(mmax,-mmin)**2/rad(i)**2&
         -w0**2/c**2*epsparr(i,nsort+2)))

   END DO

   kpermax=sqrt(maxval(kper2))
   coeffi=kpermax*(2.0_rk*a/numpt)

   IF(coeffi>khmax.OR.coeffi<khmin) THEN

      numpt=coeffi/SQRT(khmax*khmin)*numpt

      print *, ' new number of mesh points', numpt

      DEALLOCATE(rad,amatr0,amatr1p,amatr1m,bmatr,sysmatr,yy,epsparr,epsgarr,epslarr,rhsarr,p,ptmp,powm,kper2)

!-------
!Memory allocation
!-------

      ALLOCATE(rad(numpt))
      ALLOCATE(amatr1p(3*numline,numline*numpt))
      ALLOCATE(amatr1m(3*numline,numline*numpt))
      ALLOCATE(amatr0(3*numline,numline*numpt))
      ALLOCATE(bmatr(3*numline,numline*numpt))
      ALLOCATE(sysmatr(3*numline,numline*numpt))
      ALLOCATE(yy(numline*numpt))
      ALLOCATE(epsparr(numpt,nsort+2),epsgarr(numpt,nsort+2),epslarr(numpt,nsort+2))
      ALLOCATE(rhsarr(numline*numpt,mmin:mmax,nmin:nmax))
      ALLOCATE(p(numpt,nsort+2))
      ALLOCATE(ptmp(numpt,nsort+2))
      ALLOCATE(powm(numpt,nsort+2,mmin:mmax))
      ALLOCATE(kper2(numpt))

!-------
!Creating the mesh
!-------
      CALL mesh(rad)
!-------
!Preparing the right-hand-side
!-------

      DO n=nmin,nmax
         DO m=mmin,mmax

            CALL rhsmultistrap(rad,yy,m,n)

            rhsarr(:,m,n)=yy

         END DO
      END DO


!-------
!Power density array initiation
!-------
      p=0.0_rk

!-------
!Epsylon-dependent part of the matrix
!-------

      CALL projeqeps(rad,rad1,numpt1,bmatr,density,cfreq,magnf,zsort,msort,nsort)
!-------
!Matrix for power calculations
!-------
      DO i=1,numpt
         DO k=1,nsort

            CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,k,k,.FALSE.)

            epsparr(i,k)=epsp
            epsgarr(i,k)=epsg
            epslarr(i,k)=epsl

         END DO

         CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,1,0,.TRUE.)

         epsparr(i,nsort+1)=epsp
         epsgarr(i,nsort+1)=epsg
         epslarr(i,nsort+1)=epsl

         CALL epsyl(rad1,numpt1,epsp,epsg,epsl,rad(i),density,cfreq,magnf,zsort,msort,nsort,1,nsort,.TRUE.)

         epsparr(i,nsort+2)=epsp
         epsgarr(i,nsort+2)=epsg
         epslarr(i,nsort+2)=epsl

      END DO

   END IF
!-------
!Making printout for k_perp for slow and fast waves
!-------

   OPEN(unit=newunit(fidwr),file='disper.dat')										! TW20190225
!PRINT *, 'OUT disper id : ', fidwr														! TW20190225

   DO i=1,numpt
      WRITE(fidwr,'(3e15.5)') rad(i),REAL(-epslarr(i,nsort+2)/epsparr(i,nsort+2)*(mmax**2/rad(i)**2&
         -w0**2/c**2*epsparr(i,nsort+2))),&
         REAL(w0**2/c**2*epsparr(i,nsort+2)-mmax**2/rad(i)**2-nmax**2/rstar**2&
         -w0**4/c**4*epsgarr(i,nsort+2)**2/(w0**2/c**2*epsparr(i,nsort+2)-mmax**2/rad(i)**2))
   END DO
   CLOSE(fidwr)

!-------
!Start of the loop over n
!-------
   DO n=nmin,nmax

      m=0
      CALL projeq(rad,amatr0,m,n)

      m=1
      CALL projeq(rad,amatr1p,m,n)

      m=-1
      CALL projeq(rad,amatr1m,m,n)

!-------
!Making the loop over m
!-------
      CALL mcycle (n,1.0_rk,rhsarr,amatr0,amatr1p,amatr1m,bmatr,rad,epsparr,epsgarr,epslarr,nsort,p)

   END DO

!-------
!Constructing power array on the mesh for power balance corresponding to unit current
!-------


   CALL transf(rtor,rad,numpt,rad1,numpt1,nsort,p,powr)
!-------
!Calculating total power
!-------
   CALL totalpower(rad1,powr(:,nsort+2),numpt1,tp)
!-------
!Ajusting the power array
!-------
   powr=powr*(powin/tp)
   alr=2.0_rk*tp

   IF(lastcall) THEN
!-------
!Memory deallocation
!-------

      DEALLOCATE(rad,amatr0,amatr1p,amatr1m,bmatr,sysmatr,yy,epsparr,epsgarr,epslarr,rhsarr,p,ptmp,powm)

      firstcall=.TRUE.

   END IF

END SUBROUTINE RFpower
!--------------------
!--------------------
SUBROUTINE inant
!--------------------
!-------
!Reads data from the file, fills in constants in the module <antdata>.
!-------
   USE numberformat
   USE antdata
   USE rfdata
   IMPLICIT none
   INTEGER::i, fidrd, fidwr
   CHARACTER(128):: line128
   REAL(rk)::zant,dist,width1,length

   interface																									! TW20190225
      integer function newunit(unit)													! TW20190225
         integer, intent(out), optional :: unit								! TW20190225
      end function newunit														  			! TW20190225
   end interface																							! TW20190225

   OPEN (unit=newunit(fidrd),file='antparam.dat')						! TW20190225
   OPEN (unit=newunit(fidwr),file='antparamout.dat')					! TW20190225
!PRINT *, 'IN antfile id : ', fidrd												! TW20190225
!PRINT *, 'OUT antfile id : ', fidwr												! TW20190225

   READ(fidrd,*) nanten,nanthfs,sldown
   WRITE(fidwr,*) nanten,nanthfs, sldown

   READ(fidrd,*) line128
   WRITE(fidwr,*) line128

   DO i=1,nanthfs
      READ (fidrd,*) radant(i),dist,zant,width1,length,phase(i)
      WRITE (fidwr,*) radant(i),dist,zant,width1,length,phase(i)
      phiant(i)=dist/(rtor-a)
      thetaant(i)=zant/rstar
      delphiant(i)=width1/(rtor-a)
      delthetaant(i)=length/rstar

   END DO
   DO i=nanthfs+1,nanten
      READ (fidrd,*) radant(i),dist,zant,width1,length,phase(i)
      WRITE (fidwr,*) radant(i),dist,zant,width1,length,phase(i)
      phiant(i)=dist/(rtor+a)
      thetaant(i)=zant/rstar
      delphiant(i)=width1/(rtor+a)
      delthetaant(i)=length/rstar

   END DO


   CLOSE(fidrd)
   CLOSE(fidwr)

END SUBROUTINE inant
!--------------------

!--------------------
!--------------------

!--------------------
FUNCTION pointdens(x)
!--------------------
!-------
!User-defined function
!which returns the relative mesh point density as a function
!of x coordinate.
!Should be order of unity.
!-------
   USE numberformat
   USE rfdata
   IMPLICIT none
   REAL(rk)::pointdens
   REAL(rk), INTENT(IN)::x

   pointdens=1.0_rk-weight+0.56_rk*weight*a/width*exp(-((x-radac)/width)**2)

END FUNCTION pointdens
!--------------------


!--------------------
SUBROUTINE inrf
!--------------------
!-------
!Reads data from the file, fills in constants in the module <rfdata>.
!-------
   USE numberformat
   USE constants
   USE rfdata
   IMPLICIT none
   REAL(rk)::frq
   INTEGER:: fidrd,fidwr

   interface																									! TW20190225
      integer function newunit(unit)													! TW20190225
         integer, intent(out), optional :: unit								! TW20190225
      end function newunit														  			! TW20190225
   end interface																							! TW20190225

   OPEN (unit=newunit(fidrd),file='rfparamtr.dat')						! TW20190225
   OPEN (unit=newunit(fidwr),file='rfparamout.dat')					! TW20190225
!PRINT *, 'IN rffile id : ', fidrd													! TW20190225
!PRINT *, 'OUT rffile id : ', fidwr												! TW20190225

   READ(fidrd,*) numptinit
   WRITE(fidwr,*) numptinit

   READ(fidrd,*) mmin,mmax
   WRITE(fidwr,*) mmin,mmax

   READ(fidrd,*) nmin,nmax
   WRITE(fidwr,*) nmin,nmax

   READ(fidrd,*) rtor
   WRITE(fidwr,*) rtor

   READ(fidrd,*) rstar
   WRITE(fidwr,*) rstar

!READ(fidrd,*) frq
!WRITE(fidwr,*) frq
!
!w0=frq*1e6_rk*2.0_rk*pi

   READ(fidrd,*) radac
   WRITE(fidwr,*) radac

   READ(fidrd,*) width
   WRITE(fidwr,*) width

   READ(fidrd,*) weight
   WRITE(fidwr,*) weight

   READ(fidrd,*) khmax,khmin
   WRITE(fidwr,*) khmax,khmin

   READ(fidrd,*) epspmin
   WRITE(fidwr,*) epspmin

   CLOSE(fidrd)
   CLOSE(fidwr)

END SUBROUTINE inrf
!--------------------

!--------------------
SUBROUTINE mesh(rad)
!--------------------
!-------
! The subroutine constructs a 1D mesh ( array rad() )
! with the density prescribed by function pointdens().
! Besides, it marks the antenna locations
! and store these node numbers in the array antnode() placed in
! the module rfdata.
!-------
   USE numberformat
   USE rfdata
   USE antdata
   IMPLICIT none
   INTEGER,PARAMETER:: ink=10
   REAL(rk), DIMENSION(numpt), INTENT(OUT)::rad
   REAL(rk), DIMENSION(ink*numpt)::radt
   INTEGER:: i,j,jbeg,ord,count
   REAL(rk)::step,step1,coef,pointdens,rs,rleft,rright
   REAL(rk), DIMENSION(maxnumant)::rants
   REAL(rk)::lefttag,leftpre,righttag,rightpre
!-------
!-------
!Integrating point density function
!-------

   step=2.0_rk*a/(ink*numpt-1)
   radt(1)=0.0_rk

   DO i=2,ink*numpt
      rleft=rtor-a+step*(i-1)
      radt(i)=radt(i-1)+pointdens(rleft-step)/6.0_rk&
         +2.0_rk*pointdens(rleft)/3.0_rk&
         +pointdens(rleft+step)/6.0_rk

   END DO
!-------
!Scaling the result
!-------

   coef=2.0_rk*a/radt(ink*numpt)
   radt=coef*radt+rtor-a
!-------
!Forming non-uniform mesh
!-------

   step1=2.0_rk*a/(numpt-1)
   rad(1)=rtor-a
   jbeg=1

   DO i=2,numpt-1
      rs=rtor-a+step1*(i-1)
      DO j=jbeg,ink*numpt-1
         IF((radt(j)-rs)*(radt(j+1)-rs)<=0.0_rk) THEN
            rleft=rtor-a+step*(j-1)
            rright=rtor-a+step*j
            rad(i)=(rleft*(radt(j+1)-rs)+rright*(rs-radt(j)))/(radt(j+1)-radt(j))
            jbeg=j
            EXIT
         END IF
      END DO
   END DO

   rad(numpt)=rtor+a
   !-------
!Finding mesh points nearest to the straps' locations
!-------
   count=0
   DO i=2,numpt

      DO j=1,nanten
         IF((rad(i)-radant(j))*(rad(i-1)-radant(j))<=0.0_rk) THEN
            count=count+1
            IF(rad(i)-radant(j)>-rad(i-1)+radant(j)) THEN
               antnode(j)=i-1
               rants(j)=rad(i-1)
            ELSE
               antnode(j)=i
               rants(j)=rad(i)
            END IF
         END IF

      END DO
   END DO
   IF (count/=nanten) print *,' ANTNODE ERROR'

   !radant=rants
   !-------
!Ordering of straps in increasing of strap radius
!-------


   DO j=1,nanten
      ordering(j)=j
   END DO

   DO i=1,nanten
      count=0

      DO j=1,nanten-1
         IF(radant(ordering(j))>radant(ordering(j+1))) THEN
            ord=ordering(j)
            ordering(j)=ordering(j+1)
            ordering(j+1)=ord
            count=count+1
         END IF
      END DO
      IF (count==0) EXIT
   END DO

!   Open(19,file='meshpre.dat')
!
!  DO i=1,numpt
!   write(19,*) i,rad(i)
!
!  END DO
!
!  CLOSE(19)


!-------
!Puttung mesh node to first antenna location scaling 2 segments: (i) between first mesh node and antenna location and
! (ii) between antenna location and last mesh node
!-------

   lefttag=radant(ordering(1))-rad(1)
   leftpre=rad(antnode(ordering(1)))-rad(1)
   righttag=rad(numpt)-radant(ordering(1))
   rightpre=rad(numpt)-rad(antnode(ordering(1)))

   DO j=1,antnode(ordering(1))
      rad(j)=rtor-a+(rad(j)-rtor+a)*lefttag/leftpre
   END DO

   DO j=antnode(ordering(1))+1,numpt
      rad(j)=rtor+a+(rad(j)-rtor-a)*righttag/rightpre
   END DO
!-------
!Puttung mesh node to further antenna locations scaling 2 segments: (i) between previous antenna mesh node and antenna location and
! (ii) between antenna location and last mesh node
!-------


   DO i=2,nanten

      IF(antnode(ordering(i))==antnode(ordering(i-1))) CYCLE

      lefttag=radant(ordering(i))-radant(ordering(i-1))
      leftpre=rad(antnode(ordering(i)))-rad(antnode(ordering(i-1)))
      righttag=rad(numpt)-radant(ordering(i))
      rightpre=rad(numpt)-rad(antnode(ordering(i)))

      DO j=antnode(ordering(i-1)),antnode(ordering(i))
         rad(j)=radant(ordering(i-1))+(rad(j)-radant(ordering(i-1)))*lefttag/leftpre
      END DO

      DO j=antnode(ordering(i))+1,numpt
         rad(j)=rtor+a+(rad(j)-rtor-a)*righttag/rightpre
      END DO

   END DO

!     Open(19,file='meshtag.dat')
!
!  DO i=1,numpt
!   write(19,*) i,rad(i)
!
!  END DO
!
!  CLOSE(19)
!
!

END SUBROUTINE mesh


!--------------------
SUBROUTINE projeq(rad,amatr,m,n)
!--------------------
!-------
!Fills in the matrix amatr
!Uses as read only:
!rad - array of radial coordinates of grid points
!-------
   USE numberformat
   USE gaussdata
   USE rfdata
   USE eqpar
   USE constants
   IMPLICIT none
   REAL(rk), DIMENSION(numpt),INTENT(IN)::rad
   COMPLEX(rk),DIMENSION(3*numline,numline*numpt),INTENT(OUT)::amatr
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   INTEGER:: i,k,kk,kg,kv,kvv,ind1,ind2,numobs,numtag
   REAL(rk)::del,r,dla
   COMPLEX(rk)::imm,inn,ikz
   INTEGER::m,n
!-------
!indexing: arr(k,kv,kk,kvv);
!k,kk - numbers enumerating finite element functions:
!k - corresponds to row, kk - to column of the matrix.
!To learn the enumeration order see e.g. finite element
!functions <lambda12> and <lambda3>.
!kv, kvv - are introduced for vector variables being vector components:
!kv - corresponds to row, kvv - to column of the matrix.
!kv=1 - radial, kv=2 - tangent and kv=3 parallel components.
!To calculate the intermediate arrays, the Gauss quadrature formula is used.
!-------
   COMPLEX(rk), DIMENSION(nfm3,3,nfm3,3)::aa
!-------
!Functions
!-------
   REAL(rk)::lambda3,dla3dksi,lambda12


   imm=eim*m
   inn=eim*n
   ikz=inn/rstar

   amatr=(0.0_rk,0.0_rk)

!-------
!Loop over mesh segments
!-------

   DO i=1,numpt-1
      del=rad(i+1)-rad(i)

!-------
!Local matrix initiation
!-------

      aa=(0.0_rk,0.0_rk)



      DO k=1,nfm3
         DO kk=1,nfm3
            DO kg=1,kgmax
!-------
!For avoiding degeneracy
!-------

!    IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4)) THEN
               IF(i==1.AND.k==2) THEN


                  dla=lambda3(del,k,ksi(kg))

               ELSE
                  dla=lambda12(del,k,ksi(kg))
               END IF


               r=rad(i)+del*ksi(kg)

               aa(k,1,kk,1)=aa(k,1,kk,1)+weigh(kg)*dla &
                  *lambda3(del,kk,ksi(kg))*(-ikz**2*r+(m**2)/r)*del

               aa(k,1,kk,2)=aa(k,1,kk,2)+weigh(kg)*dla &
                  *imm*(lambda3(del,kk,ksi(kg))/r+dla3dksi(del,kk,ksi(kg))/del)*del

               aa(k,1,kk,3)=aa(k,1,kk,3)+weigh(kg)*dla &
                  *ikz*r*dla3dksi(del,kk,ksi(kg))

               aa(k,2,kk,1)=aa(k,2,kk,1)-weigh(kg)*lambda3(del,kk,ksi(kg))&
                  *imm*dla3dksi(del,k,ksi(kg))/r

               aa(k,2,kk,2)=aa(k,2,kk,2)+weigh(kg) &
                  *(-lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del &
                  *ikz**2+lambda3(del,kk,ksi(kg))*dla3dksi(del,k,ksi(kg))/r &
                  +dla3dksi(del,kk,ksi(kg))*dla3dksi(del,k,ksi(kg))/del)


               aa(k,2,kk,3)=aa(k,2,kk,3)+weigh(kg)*imm*ikz &
                  *lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del/r

               aa(k,3,kk,1)=aa(k,3,kk,1)-weigh(kg)*lambda3(del,kk,ksi(kg))&
                  *ikz*r*dla3dksi(del,k,ksi(kg))

               aa(k,3,kk,2)=aa(k,3,kk,2)+weigh(kg)*imm*ikz &
                  *lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del

               aa(k,3,kk,3)=aa(k,3,kk,3)+weigh(kg)*(m**2/r &
                  *lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del &
                  +dla3dksi(del,kk,ksi(kg))*dla3dksi(del,k,ksi(kg))*r/del)

            END DO
         END DO
      END DO
!-------
!filling up global matrix
!-------
      numobs=1
      numtag=1

      DO k=1,nfm3
         DO kk=1,nfm3
            DO kv=1,3
               DO kvv=1,3

                  CALL indexes(k,kv,kk,kvv,numobs,numtag,i,ind1,ind2)
                  amatr(ind1,ind2)=amatr(ind1,ind2)+aa(k,kv,kk,kvv)
               END DO
            END DO
         END DO
      END DO

   END DO

END SUBROUTINE projeq
!--------------------


!--------------------
SUBROUTINE projeqeps(rad,rad1,numpt1,bmatr,density,cfreq,magnf,zsort,msort,nsort)
!--------------------
!-------
!Fills in the matrix bmatr
!Uses as read only:
!rad - array of radial coordinates of grid points;
!-------
   USE numberformat
   USE gaussdata
   USE rfdata
   USE eqpar
   USE constants
   IMPLICIT none
   REAL(rk), DIMENSION(numpt),INTENT(IN)::rad
   COMPLEX(rk),DIMENSION(3*numline,numline*numpt),INTENT(OUT)::bmatr
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   INTEGER:: i,k,kk,kg,kv,kvv,ind1,ind2,numobs,numtag
   INTEGER,INTENT(IN):: nsort
   REAL(rk)::del,r,dla,d2la
   COMPLEX(rk)::epsp,epsg,epsl
   INTEGER,INTENT(IN)::numpt1
   REAL(rk),DIMENSION(2*numpt1,nsort)::density
   REAL(rk),DIMENSION(numpt1,nsort+1)::cfreq
   REAL(rk),DIMENSION(numpt1)::rad1
   REAL(rk),DIMENSION(numpt1)::magnf
   INTEGER,DIMENSION(nsort)::zsort,msort
!-------
!indexing: arr(k,kv,kk,kvv);
!k,kk - numbers enumerating finite element functions:
!k - corresponds to row, kk - to column of the matrix.
!To learn the enumeration order see e.g. finite element
!functions <lambda2> and <lambda3>.
!kv, kvv - are introduced for vector variables being vector components:
!kv - corresponds to row, kvv - to column of the matrix.
!kv=1 - radial, kv=2 - tangent and kv=3 parallel components.
!To calculate the intermediate arrays, the Gauss quadrature formula is used.
!-------
   COMPLEX(rk), DIMENSION(nfm3,3,nfm3,3)::aa
!-------
!Functions
!-------
   REAL(rk)::lambda3,dla3dksi,lambda12,dla12dksi


   bmatr=(0.0_rk,0.0_rk)

!-------
!Loop over mesh segments
!-------

   DO i=1,numpt-1
      del=rad(i+1)-rad(i)


      aa=(0.0_rk,0.0_rk)

      DO k=1,nfm3
         DO kk=1,nfm3
            DO kg=1,kgmax
               r=rad(i)+del*ksi(kg)

               !IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4)) THEN
               IF(i==1.AND.k==2) THEN

                  dla=lambda3(del,k,ksi(kg))
                  d2la=dla3dksi(del,k,ksi(kg))
               ELSE
                  dla=lambda12(del,k,ksi(kg))
                  d2la=dla12dksi(del,k,ksi(kg))
               END IF

               CALL  epsyl(rad1,numpt1,epsp,epsg,epsl,r,density,cfreq,magnf,zsort,msort,nsort,1,nsort,.TRUE.)


               aa(k,1,kk,1)=aa(k,1,kk,1)+weigh(kg)*dla &
                  *lambda3(del,kk,ksi(kg))*epsp*r/c**2*del

               aa(k,1,kk,3)=aa(k,1,kk,3)-weigh(kg)*dla &
                  *lambda3(del,kk,ksi(kg))*eim*epsg*r/c**2*del

               aa(k,3,kk,1)=aa(k,3,kk,1)+weigh(kg)*lambda3(del,kk,ksi(kg))&
                  *lambda3(del,k,ksi(kg))*eim*epsg*del*r/c**2

               aa(k,3,kk,3)=aa(k,3,kk,3)+weigh(kg) &
                  *lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del &
                  *epsp*r/c**2

               aa(k,2,kk,2)=aa(k,2,kk,2)+weigh(kg) &
                  *lambda3(del,k,ksi(kg))*lambda3(del,kk,ksi(kg))*del &
                  *epsl/c**2



            END DO
         END DO
      END DO


      numobs=1
      numtag=1

      DO k=1,nfm3
         DO kk=1,nfm3
            DO kv=1,3
               DO kvv=1,3

                  CALL indexes(k,kv,kk,kvv,numobs,numtag,i,ind1,ind2)
                  bmatr(ind1,ind2)=bmatr(ind1,ind2)+aa(k,kv,kk,kvv)
               END DO
            END DO
         END DO
      END DO



   END DO

END SUBROUTINE projeqeps



!--------------------
SUBROUTINE singleharm (m,n,curamp,rhsarr,amatr0,amatr1p,amatr1m,bmatr,rad,epsparr,epsgarr,epslarr,nsort,powz)
!--------------------
! Calculates power deposition profile for single Fourier harmonic
!-------
   USE numberformat
   USE eqpar
   USE rfdata
   IMPLICIT none
   INTEGER,INTENT(IN)::m,n,nsort
   REAL(rk)::curamp
   REAL(rk), DIMENSION(numpt)::rad
   COMPLEX(rk),DIMENSION(numline*numpt,mmin:mmax,nmin:nmax)::rhsarr
   COMPLEX(rk),DIMENSION(numpt,nsort+2)::epsparr,epsgarr,epslarr
   COMPLEX(rk),DIMENSION(3*numline,numline*numpt)::sysmatr,amatr1p,amatr1m,amatr0,bmatr
   COMPLEX(rk),DIMENSION(numline*numpt)::yy
   REAL(rk),DIMENSION(numpt,nsort+2)::powz

   INTEGER::i


   yy=curamp*rhsarr(:,m,n)


   sysmatr=(1.0_rk-m**2)*amatr0+(m+m**2)/2.0_rk*amatr1p+(-m+m**2)/2.0_rk*amatr1m-w0**2*bmatr

   CALL bndcon(rad,sysmatr,epsparr,epsgarr,epslarr,nsort,n)


   CALL clubl3diag(sysmatr,numline,numpt)

   CALL cslvbl3diag(sysmatr,yy,numline,numpt)

   CALL power(epsparr,epsgarr,epslarr,nsort,rad,yy,powz)

!OPEN (33,file='fields.dat')
!DO i=1,numpt
!
! WRITE(33,'(8e15.5)') rad(i), yy(numline*(i-1)+1),yy(numline*(i-1)+numfem+1),yy(numline*(i-1)+2*numfem+1)
!
!END DO
!CLOSE(33)
!

END SUBROUTINE singleharm
!--------------------


!--------------------
SUBROUTINE mcycle (n,curamp,rhsarr,amatr0,amatr1p,amatr1m,bmatr,rad,epsparr,epsgarr,epslarr,nsort,powz)
!--------------------
!-------
!Loop over m with OpenMP directives
!-------
!-------
   USE numberformat
   USE eqpar
   USE rfdata
   IMPLICIT none
   INTEGER,INTENT(IN)::n,nsort
   REAL(rk)::curamp
   REAL(rk), DIMENSION(numpt)::rad
   COMPLEX(rk),DIMENSION(numline*numpt,mmin:mmax,nmin:nmax)::rhsarr
   COMPLEX(rk),DIMENSION(numpt,nsort+2)::epsparr,epsgarr,epslarr
   COMPLEX(rk),DIMENSION(3*numline,numline*numpt)::amatr1p,amatr1m,amatr0,bmatr
   REAL(rk),DIMENSION(numpt,nsort+2)::powz
   REAL(rk),DIMENSION(numpt,nsort+2)::powa,powb
   INTEGER::m

   powb=0.0_rk
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(m,powa),REDUCTION(+:powb)


   DO m=mmin,mmax

      CALL singleharm (m,n,curamp,rhsarr,amatr0,amatr1p,amatr1m,bmatr,rad,epsparr,epsgarr,epslarr,nsort,powa)
! CALL dummm (m,n,nsort,powa)

      powb=powb+powa

   END DO

!$OMP END PARALLEL DO

   powz=powz+powb

END SUBROUTINE mcycle
!--------------------


!--------------------
SUBROUTINE bndcon(rad,amatr,epsparr,epsgarr,epslarr,nsort,n)
!--------------------
!-------
!Imposes the boundary conditions for matrix
!-------
   USE numberformat
   USE eqpar
   USE rfdata
   IMPLICIT none

   REAL(rk),INTENT(IN):: rad(numpt)
   COMPLEX(rk),DIMENSION(3*numline,numline*numpt),INTENT(INOUT)::amatr
   COMPLEX(rk),DIMENSION(numpt,nsort+2)::epsparr,epsgarr,epslarr

   INTEGER::kf,n,nsort

!-------
!Boundary conditions nullifying E_tau and E_// at the boundary
!-------
   kf=numline*(numpt-1)
   amatr(:,kf+numfem+1)=(0.0_rk,0.0_rk)
   amatr(:,kf+2*numfem+1)=(0.0_rk,0.0_rk)

   amatr(numline+numfem+1,kf+numfem+1)=(1.0_rk,0.0_rk)
   amatr(numline+2*numfem+1,kf+2*numfem+1)=(1.0_rk,0.0_rk)

!amatr(:,kf+2)=(0.0_rk,0.0_rk)
!amatr(numline+2,kf+2)=epsparr(numpt,nsort+2)
!amatr(numline+1,kf+2)=(epsparr(numpt,nsort+2)*rad(numpt)-epsparr(numpt-1,nsort+2)*rad(numpt-1))/(rad(numpt)-rad(numpt-1))/rad(numpt)&
!                      +n*epsgarr(numpt,nsort+2)/rstar
!amatr(numline+4,kf+2)=(0.0_rk,1.0_rk)*epsgarr(numpt,nsort+2)


   kf=0
   amatr(:,kf+numfem+1)=(0.0_rk,0.0_rk)
   amatr(:,kf+2*numfem+1)=(0.0_rk,0.0_rk)

   amatr(numline+numfem+1,kf+numfem+1)=(1.0_rk,0.0_rk)
   amatr(numline+2*numfem+1,kf+2*numfem+1)=(1.0_rk,0.0_rk)

!amatr(:,kf+2)=(0.0_rk,0.0_rk)
!amatr(numline+2,kf+2)=epsparr(1,nsort+2)
!amatr(numline+1,kf+2)=(epsparr(2,nsort+2)*rad(2)-epsparr(1,nsort+2)*rad(1))/(rad(2)-rad(1))/rad(1)&
!+n*epsgarr(1,nsort+2)/rstar
!amatr(numline+4,kf+2)=(0.0_rk,1.0_rk)*epsgarr(1,nsort+2)


END SUBROUTINE bndcon
!--------------------



!--------------------
FUNCTION lambda3(del,k,ksi)
!--------------------
!-------
! Calculates the value of hat function at the interval
!-------
   USE numberformat
   IMPLICIT none
   REAL(rk), INTENT(IN):: del,ksi
   REAL(rk):: lambda3
   INTEGER, INTENT(IN):: k

   SELECT CASE(k)

    CASE(1)
      lambda3=(1.0_rk-ksi)**2*(1.0_rk+2.0_rk*ksi)

    CASE(2)
      lambda3=del*(1.0_rk-ksi)**2*ksi

    CASE(3)
      lambda3=ksi**2*(3.0_rk-2.0_rk*ksi)

    CASE(4)
      lambda3=-del*(1.0_rk-ksi)*ksi**2

    CASE DEFAULT

      STOP ' Lambda3: hat function number is wrong'

   END SELECT

END FUNCTION lambda3
!--------------------



!--------------------
FUNCTION dla3dksi(del,k,ksi)
!--------------------
!-------
! Calculates the derivatinve of hat function at the interval
!-------
   USE numberformat
   IMPLICIT none
   REAL(rk), INTENT(IN):: del,ksi
   REAL(rk):: dla3dksi
   INTEGER, INTENT(IN):: k

   SELECT CASE(k)

    CASE(1)
      dla3dksi=-6.0_rk*ksi*(1.0_rk-ksi)

    CASE(2)
      dla3dksi=del*(1.0_rk-ksi)*(1.0_rk-3.0_rk*ksi)

    CASE(3)
      dla3dksi=6.0_rk*ksi*(1.0_rk-ksi)

    CASE(4)
      dla3dksi=del*ksi*(3.0_rk*ksi-2.0_rk)

    CASE DEFAULT

      STOP ' dla3dksi: hat function number is wrong'

   END SELECT

END FUNCTION dla3dksi
!--------------------


!--------------------
FUNCTION lambda12(del,k,ksi)
!--------------------
!-------
! Calculates the value of hat function at the interval
!-------
   USE numberformat
   IMPLICIT none
   REAL(rk), INTENT(IN):: del,ksi
   REAL(rk):: lambda12
   INTEGER, INTENT(IN):: k

   SELECT CASE(k)

    CASE(1)
      lambda12=1.0_rk-ksi

    CASE(2)
      lambda12=del*(1.0_rk-ksi)*ksi

    CASE(3)
      lambda12=ksi

    CASE(4)
      lambda12=-del*(1.0_rk-ksi)*ksi

    CASE DEFAULT

      STOP ' Lambda12: hat function number is wrong'

   END SELECT

END FUNCTION lambda12
!--------------------



!--------------------
FUNCTION dla12dksi(del,k,ksi)
!--------------------
!-------
! Calculates the value of hat function at the interval
!-------
   USE numberformat
   IMPLICIT none
   REAL(rk), INTENT(IN):: del,ksi
   REAL(rk):: dla12dksi
   INTEGER, INTENT(IN):: k

   SELECT CASE(k)

    CASE(1)
      dla12dksi=-1.0_rk

    CASE(2)
      dla12dksi=del*(1.0_rk-2.0_rk*ksi)

    CASE(3)
      dla12dksi=1.0_rk

    CASE(4)
      dla12dksi=-del*(1.0_rk-2.0_rk*ksi)

    CASE DEFAULT

      STOP ' dla12dksi: hat function number is wrong'

   END SELECT

END FUNCTION dla12dksi
!--------------------



!--------------------
SUBROUTINE indexes(k,kv,kk,kvv,neqobs,neqlink,i,ind1,ind2)
!--------------------
!-------
! Calculates the  indexes of matrix
!
!k,kk - numbers enumerating finite element functions:
!k - corresponds to row, kk - to column of the matrix.
!To learn the enumeration order see e.g. finite element
!functions <lambda2> and <lambda3>.
!kv, kvv - are  vector components:
!kv - corresponds to row, kvv - to column of the matrix.
!kv=1 - radial, kv=2 - tangent and kv=3 parallel components.
! For scalar variables kv=1 or kvv=1.
!neqobs - the number of equation,neqlink  - the number of variable
!which contributes to the equation. For vector variables
!<neqobs> or <neqlink> are the numbers of the radial components.
!i - number of the segment.
!ind1,ind2 -first and second index of matrix array <sysmatr>
!The composition of the vector of the system is the following:
!...
!<finite element 1, variable n, node i>
!<finite element 2, variable n, node i>
!<finite element 1, variable n+1, node i>
!...
!Matrix of the system a(i,j) is of size (numline*numpt)x(numline*numpt);
!it is stored in array sysmatr(k1,k2),so that k2=i. k1=j-Funct(i).
!
!Funct(i)=-numline for i=1...numline;
!Funct(i)=0, for i=numline+1...2*numline,
!... Funct(i)=(s-1)*numline, for i=s*numlin+1...(s+1)*numline.
!-------
   USE numberformat
   USE eqpar
   IMPLICIT none
   INTEGER,INTENT(IN)::k,kv,kk,kvv,neqobs,neqlink,i
   INTEGER,INTENT(OUT)::ind1,ind2

   ind2=numline*(i-1+(k-1)/numfem)+numfem*(neqobs+kv-2) &
      +k-numfem*((k-1)/numfem)

   ind1=numline*(i-1+(kk-1)/numfem)+numfem*(neqlink+kvv-2) &
      +kk-numfem*((kk-1)/numfem) &
      +numline-numline*((ind2-1)/numline)

END SUBROUTINE indexes
!--------------------



!--------------------
SUBROUTINE clubl3diag(a,nb,nl)
!--------------------
!-------
!performs LU decomposition L*U=a for complex block-3diagonal matrix <a>
!the result  is stored in <a>
!nb - size of the block
!nl - number of block equations
!-------
   USE numberformat
   IMPLICIT none
   INTEGER, INTENT(IN)::nb,nl
   COMPLEX(rk),INTENT(INOUT):: a(3*nb,nb*nl)
   COMPLEX(rk):: sl,sc
   INTEGER::i,j,k,ii,indeks,indi,indj

   indeks(ii)=nb*((ii-1)/nb)-nb

   DO i=1,nb*nl
      indi=indeks(i)
      DO j=max(indi+1,1),i-1
         indj=indeks(j)
         sc=a(i-indj,j)

         DO k=max(indi+1,1),j-1
            sc=sc-a(k-indj,j)*a(i-indeks(k),k)
         END DO
         a(i-indj,j)=sc/a(j-indj,j)
      END DO


      DO j=max(indi+1,1),i
         indj=indeks(j)
         sl=a(j-indi,i)

         DO k=max(indi+1,1),j-1
            sl=sl-a(k-indi,i)*a(j-indeks(k),k)
         END DO

         a(j-indi,i)=sl
      END DO

   END DO

END SUBROUTINE clubl3diag
!--------------------



!--------------------
SUBROUTINE cslvbl3diag(a,y,nb,nl)
!--------------------
!-------
!solves linear equations system a*x=y for complex block-3diagonal matrix <a>
!the result <x> is stored in <y>
!nb - size of the block
!nl - number of block equations
!-------
   USE numberformat
   IMPLICIT none
   INTEGER, INTENT(IN)::nb,nl
   COMPLEX(rk),INTENT(IN):: a(3*nb,nb*nl)
   COMPLEX(rk),INTENT(INOUT):: y(nb*nl)
   COMPLEX(rk):: s
   INTEGER::i,j,ii,indeks,indi

   indeks(ii)=nb*((ii-1)/nb)-nb

   DO i=1,nb*nl
      indi=indeks(i)
      s=(0.0_rk,0.0_rk)
      DO j=max(indi+1,1),i-1
         s=s+a(j-indi,i)*y(j)
      END DO
      y(i)=(y(i)-s)/a(i-indi,i)

   END DO


   DO i=nb*nl,1,-1
      indi=indeks(i)
      s=(0.0_rk,0.0_rk)
      DO j=min(indi+3*nb,nb*nl),i+1,-1
         s=s+a(j-indi,i)*y(j)
      END DO
      y(i)=y(i)-s

   END DO

END SUBROUTINE cslvbl3diag
!--------------------


!----------------
SUBROUTINE rhsmultistrap(rad,yy,m,n)
!--------------------
!------
!Calculation the right hand side of Maxwell equation
!Calculating the antenna currents
!i*w0*mu0*j
!------
!Four strap antenna
   USE numberformat
   USE constants
   USE gaussdata
   USE rfdata
   USE eqpar
   USE antdata
   IMPLICIT NONE
!INTEGER, PARAMETER::numpt=101
   REAL(rk),INTENT(IN)::rad(numpt)
   COMPLEX(rk),DIMENSION(numline*numpt),INTENT(OUT)::yy
   INTEGER,INTENT(IN)::m,n
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   INTEGER:: i,k,kg,j
   REAL(rk)::lambda3,lambda12
!REAL(rk)::r,del,dla,pi,mu0,totcur
   REAL(rk)::r,del,dla,totcur,k0rstar
   INTEGER::kf
   totcur=1.0_rk
   k0rstar=sldown*w0/c*rstar



   yy=(0.0_rk,0.0_rk)

   DO j=1,nanthfs
      DO i=1,antnode(j)-1
         del=rad(i+1)-rad(i)
         DO k=1,nfm3
            DO kg=1,kgmax


               r=rad(i)+del*ksi(kg)
               ! IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4)) THEN
               IF(i==1.AND.k==2) THEN
                  dla=lambda3(del,k,ksi(kg))
               ELSE
                  dla=lambda12(del,k,ksi(kg))
               END IF


!left node
               IF(k<=2) THEN

                  IF(m/=0) THEN

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF

               END IF
!right node
               IF(k>2) THEN

                  IF(m/=0) THEN

                     yy(numline*i+k-2)=yy(numline*i+k-2)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*i+k-2)=yy(numline*i+k-2)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF
               END IF


            END DO
         END DO
      END DO

      i=antnode(j)
      IF(m==0) THEN

         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/(2.0_rk*pi**2)*totcur/&
            delphiant(j)/cos(k0rstar*delthetaant(j)/2.0_rk)&
            *exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            (sin((n-k0rstar)*delthetaant(j)/2.0_rk)/(n-k0rstar)+sin((n+k0rstar)*delthetaant(j)/2.0_rk)/(n+k0rstar))&
            *delphiant(j)/2.0_rk



      ELSE
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/(2.0_rk*pi**2)*totcur/&
            delphiant(j)/cos(k0rstar*delthetaant(j)/2.0_rk)&
            *exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            (sin((n-k0rstar)*delthetaant(j)/2.0_rk)/(n-k0rstar)+sin((n+k0rstar)*delthetaant(j)/2.0_rk)/(n+k0rstar))&
            *sin(m*delphiant(j)/2.0_rk)/m
      END IF



   END DO

   DO j=nanthfs+1,nanten

      DO i=antnode(j),numpt-1
         del=rad(i+1)-rad(i)
         DO k=1,nfm3
            DO kg=1,kgmax


               r=rad(i)+del*ksi(kg)
               !IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4))
               IF(i==1.AND.k==2) THEN


                  dla=lambda3(del,k,ksi(kg))
               ELSE
                  dla=lambda12(del,k,ksi(kg))
               END IF


!left node
               IF(k<=2) THEN

                  IF(m/=0) THEN

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF

               END IF
!right node
               IF(k>2) THEN

                  IF(m/=0) THEN

                     yy(numline*i+k-2)=yy(numline*i+k-2)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*i+k-2)=yy(numline*i+k-2)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF
               END IF


            END DO
         END DO
      END DO

      i=antnode(j)
      IF(m==0) THEN

         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/(2.0_rk*pi**2)*totcur/&
            delphiant(j)/cos(k0rstar*delthetaant(j)/2.0_rk)&
            *exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            (sin((n-k0rstar)*delthetaant(j)/2.0_rk)/(n-k0rstar)+sin((n+k0rstar)*delthetaant(j)/2.0_rk)/(n+k0rstar))&
            *delphiant(j)/2.0_rk



      ELSE
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/(2.0_rk*pi**2)*totcur/&
            delphiant(j)/cos(k0rstar*delthetaant(j)/2.0_rk)&
            *exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            (sin((n-k0rstar)*delthetaant(j)/2.0_rk)/(n-k0rstar)+sin((n+k0rstar)*delthetaant(j)/2.0_rk)/(n+k0rstar))&
            *sin(m*delphiant(j)/2.0_rk)/m
      END IF



   END DO

   kf=0
!yy(kf+2)=(0.0_rk,0.0_rk)
   yy(kf+numfem+1)=(0.0_rk,0.0_rk)
   yy(kf+2*numfem+1)=(0.0_rk,0.0_rk)

   kf=numline*(numpt-1)

!yy(kf+2)=(0.0_rk,0.0_rk)
   yy(kf+numfem+1)=(0.0_rk,0.0_rk)
   yy(kf+2*numfem+1)=(0.0_rk,0.0_rk)

END SUBROUTINE rhsmultistrap
!-------------------

!----------------
SUBROUTINE rhsmultistrap1(rad,yy,m,n)
!--------------------
!------
!Calculation the right hand side of Maxwell equation
!Calculating the antenna currents
!i*w0*mu0*j
!------
!Four strap antenna
   USE numberformat
   USE constants
   USE gaussdata
   USE rfdata
   USE eqpar
   USE antdata
   IMPLICIT NONE
!INTEGER, PARAMETER::numpt=101
   REAL(rk),INTENT(IN)::rad(numpt)
   COMPLEX(rk),DIMENSION(numline*numpt),INTENT(OUT)::yy
   INTEGER,INTENT(IN)::m,n
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   INTEGER:: i,k,kg,j
   REAL(rk)::lambda3,lambda12
!REAL(rk)::r,del,dla,pi,mu0,totcur
   REAL(rk)::r,del,dla,totcur
   INTEGER::kf
! I=1 -total current
!pi=3.14154_rk
!rtor=155.0_rk
!rstar=300.0_rk
!mu0=4.0_rk*pi*1e-9_rk
   totcur=1.0_rk

   yy=(0.0_rk,0.0_rk)

   DO j=1,nanthfs
      DO i=1,antnode(j)-1
         del=rad(i+1)-rad(i)
         DO k=1,nfm3
            DO kg=1,kgmax


               r=rad(i)+del*ksi(kg)
               !IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4)) THEN
               IF(i==1.AND.k==2) THEN

                  dla=lambda3(del,k,ksi(kg))
               ELSE
                  dla=lambda12(del,k,ksi(kg))
               END IF


!left node
               IF(k<=2) THEN

                  IF(m/=0) THEN

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF

               END IF
!right node
               IF(k>2) THEN

                  IF(m/=0) THEN

                     yy(numline*i+k-2)=yy(numline*i+k-2)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*i+k-2)=yy(numline*i+k-2)-&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF
               END IF


            END DO
         END DO
      END DO

      i=antnode(j)
      IF(m==0.AND.n==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            delthetaant(j)/2.0_rk*delphiant(j)/2.0_rk
      ELSE IF(m==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            sin(n*delthetaant(j)/2.0_rk)/n*delphiant(j)/2.0_rk


      ELSE IF(n==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            delthetaant(j)/2.0_rk*sin(m*delphiant(j)/2.0_rk)/m


      ELSE
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            sin(n*delthetaant(j)/2.0_rk)/n*sin(m*delphiant(j)/2.0_rk)/m
      END IF



   END DO

   DO j=nanthfs+1,nanten

      DO i=antnode(j),numpt-1
         del=rad(i+1)-rad(i)
         DO k=1,nfm3
            DO kg=1,kgmax


               r=rad(i)+del*ksi(kg)
               !IF((i==1.AND.k==2).OR.(i==numpt-1.AND.k==4)) THEN
               IF(i==1.AND.k==2) THEN

                  dla=lambda3(del,k,ksi(kg))
               ELSE
                  dla=lambda12(del,k,ksi(kg))
               END IF


!left node
               IF(k<=2) THEN

                  IF(m/=0) THEN

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*(i-1)+k)=yy(numline*(i-1)+k)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF

               END IF
!right node
               IF(k>2) THEN

                  IF(m/=0) THEN

                     yy(numline*i+k-2)=yy(numline*i+k-2)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*sin(m*delphiant(j)/2.0_rk)/m
                  ELSE

                     yy(numline*i+k-2)=yy(numline*i+k-2)+&
                        phase(j)*w0*mu0*weigh(kg)*dla*del/pi**2*totcur/&
                        delphiant(j)/rstar*exp(-eim*n*thetaant(j))*&
                        sin(n*delthetaant(j)/2.0_rk)*delphiant(j)/2.0_rk
                  END IF
               END IF


            END DO
         END DO
      END DO

      i=antnode(j)
      IF(m==0.AND.n==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            delthetaant(j)/2.0_rk*delphiant(j)/2.0_rk
      ELSE IF(m==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            sin(n*delthetaant(j)/2.0_rk)/n*delphiant(j)/2.0_rk


      ELSE IF(n==0) THEN
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            delthetaant(j)/2.0_rk*sin(m*delphiant(j)/2.0_rk)/m


      ELSE
         yy(numline*(i-1)+2*numfem+1)=yy(numline*(i-1)+2*numfem+1)+&
            phase(j)*eim*w0*mu0/pi**2*totcur/&
            delphiant(j)*exp(-eim*n*thetaant(j)-eim*m*phiant(j))*&
            sin(n*delthetaant(j)/2.0_rk)/n*sin(m*delphiant(j)/2.0_rk)/m
      END IF



   END DO

   kf=0
   yy(kf+numfem+1)=(0.0_rk,0.0_rk)
   yy(kf+2*numfem+1)=(0.0_rk,0.0_rk)

   kf=numline*(numpt-1)

   yy(kf+numfem+1)=(0.0_rk,0.0_rk)
   yy(kf+2*numfem+1)=(0.0_rk,0.0_rk)

END SUBROUTINE rhsmultistrap1
!-------------------

!-------------------
SUBROUTINE power(epsparr,epsgarr,epslarr,nsort,rad,yy,pows)
!--------------------
!-------
!Calculating the power density
!-------
   USE numberformat
   USE gaussdata
   USE rfdata
   USE constants
   USE eqpar
   IMPLICIT none
   INTEGER,INTENT(IN)::nsort
   COMPLEX(rk),DIMENSION(numpt,nsort+2),INTENT(IN)::epsparr,epsgarr,epslarr
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   REAL(rk),INTENT(IN)::rad(numpt)
   COMPLEX(rk),DIMENSION(numline*numpt),INTENT(IN)::yy
   REAL(rk),DIMENSION(numpt,nsort+2)::pows
   INTEGER::i,k

   pows=0.0_rk

! DO i=1,numpt-1
   DO i=1,numpt
      DO k=1,nsort

         pows(i,k)=REAL(w0)*eps0/2.0_rk*IMAG(ABS(yy(numline*(i-1)+1))**2*epsparr(i,k)+&
            ABS(yy(numline*(i-1)+2*numfem+1))**2*epsparr(i,k)+&
            ABS(yy(numline*(i-1)+numfem+1))**2*epslarr(i,k)-&
            2.0_rk*yy(numline*(i-1)+1)*IMAG(epsgarr(i,k))*&
            CONJG(yy(numline*(i-1)+2*numfem+1)))

      END DO

      pows(i,nsort+1)=REAL(w0)*eps0/2.0_rk*IMAG(ABS(yy(numline*(i-1)+1))**2*epsparr(i,nsort+1)+&
         ABS(yy(numline*(i-1)+2*numfem+1))**2*epsparr(i,nsort+1)+&
         ABS(yy(numline*(i-1)+numfem+1))**2*epslarr(i,nsort+1)-&
         2.0_rk*yy(numline*(i-1)+1)*IMAG(epsgarr(i,nsort+1))*&
         CONJG(yy(numline*(i-1)+2*numfem+1)))

      pows(i,nsort+2)=REAL(w0)*eps0/2.0_rk*IMAG(ABS(yy(numline*(i-1)+1))**2*epsparr(i,nsort+2)+&
         ABS(yy(numline*(i-1)+2*numfem+1))**2*epsparr(i,nsort+2)+&
         ABS(yy(numline*(i-1)+numfem+1))**2*epslarr(i,nsort+2)-&
         2.0_rk*yy(numline*(i-1)+1)*IMAG(epsgarr(i,nsort+2))*&
         CONJG(yy(numline*(i-1)+2*numfem+1)))

   END DO

END SUBROUTINE power
!--------------------



!--------------------
SUBROUTINE epsyl(rad1,numpt1,epsp,epsg,epsl,r,density,cfreq,magnf,zsort,msort,nsort,nsortst,nsortfin,truel)
!--------------------
! Calculates dielectric tensor components at radius r
!-------
   USE numberformat
   USE constants
   USE rfdata
   IMPLICIT none
   REAL(rk)::corr
   REAL(rk)::r
   INTEGER,INTENT(IN)::numpt1
   REAL(rk),DIMENSION(2*numpt1,nsort)::density
   REAL(rk),DIMENSION(numpt1,nsort+1)::cfreq
   INTEGER,INTENT(IN)::nsort
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::rad1
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::magnf
   COMPLEX(rk), PARAMETER::eim=(0.0_rk,1.0_rk)
   COMPLEX(rk)::ompi2,ompe2,omci,omce,freqin,freqen,freqei,epsl,epsp,epsg
   INTEGER,INTENT(IN)::nsortst,nsortfin
   LOGICAL,INTENT(IN)::truel
   INTEGER,DIMENSION(nsort)::zsort,msort
   INTEGER::i
   REAL(rk)::bval,dense
   REAL(rk),DIMENSION(nsort)::densval
   REAL(rk),DIMENSION(nsort+1)::freqval

   epsp=1.0_rk

   epsg=0.0_rk

   epsl=1.0_rk

   CALL values(rad1,numpt1,r,density,cfreq,nsort,magnf,densval,freqval,bval)
!--------------------
! Ion contributions
!-------

   DO i=nsortst,nsortfin

      ompi2=4.0_rk*pi*ec**2*zsort(i)**2*densval(i)/pm/msort(i)

      omci=ec*bval*zsort(i)/(pm*msort(i)*c)


      freqin=eim*freqval(i)


      epsp=epsp+ompi2/(omci**2-(w0+freqin)**2)*(w0+freqin)/w0

      epsg=epsg+ompi2/(omci**2-(w0+freqin)**2)*omci/w0

   END DO
!--------------------
! Electron contribution
!-------

   IF (truel) THEN

      dense=0.0_rk

      DO i=1,nsort

         dense=dense+densval(i)*zsort(i)

      END DO

      freqei=0.0_rk

      ompe2=4.0_rk*pi*ec**2*dense/em

      omce=ec*bval/(em*c)

      freqen=eim*freqval(nsort+1)


      epsp=epsp+ompe2/omce**2*(w0+freqen)/w0

      epsg=epsg-ompe2/omce/w0

      epsl=epsl-ompe2/(w0*(w0+freqen+freqei))

   END IF

!--------------------
! Penalty method  contribution
!-------
   IF (nsortst==1.AND.nsortfin==nsort.AND.truel) THEN

      !epspmin=0.2_rk

      IF(ABS(epsp)<epspmin) THEN

         corr=SQRT(epspmin**2-REAL(epsp)**2)-IMAG(epsp)

         !corr=SQRT((REAL(epsp)**2/2.0_rk/epspmin+epspmin/2.0_rk)**2-REAL(epsp)**2)-IMAG(epsp)

         corr=SQRT(IMAG(epsp)**2-ABS(epsp)**2/2.0_rk+epspmin**2/2.0_rk)-IMAG(epsp)

         epsp=epsp+CMPLX(0.0_rk,corr)

      END IF

   END IF

END SUBROUTINE epsyl
!--------------------


!--------------------
SUBROUTINE values(rad1,numpt1,r,density,cfreq,nsort,magnf,densval,freqval,bval)
!--------------------
! Calculates values at point r using interpolation
!-------
   USE numberformat
   USE rfdata
   USE constants
   IMPLICIT none
   REAL(rk)::r
   INTEGER,INTENT(IN)::nsort
   REAL(rk),DIMENSION(nsort),INTENT(OUT)::densval
   REAL(rk),DIMENSION(nsort+1),INTENT(OUT)::freqval
   REAL(rk),INTENT(OUT)::bval
   INTEGER::i1
   INTEGER,INTENT(IN)::numpt1
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::rad1
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::magnf
   REAL(rk),DIMENSION(2*numpt1,nsort),INTENT(IN)::density
   REAL(rk),DIMENSION(numpt1,nsort+1),INTENT(IN)::cfreq
   REAL(rk):: ksi,del,lambda3

   DO i1=1,numpt1-1

      IF (r>=rad1(i1)+rtor.and.r<=rad1(i1+1)+rtor) EXIT

   END DO

   IF (i1==numpt1) print *, 'values:point out range'

   del=rad1(i1+1)-rad1(i1)
   ksi=(r-rad1(i1)-rtor)/del

   densval(:)=density(2*i1-1,:)*lambda3(del,1,ksi)+&
      density(2*i1,:)*lambda3(del,2,ksi)+&
      density(2*i1+1,:)*lambda3(del,3,ksi)+&
      density(2*i1+2,:)*lambda3(del,4,ksi)

   freqval(:)=cfreq(i1,:)*(1.0_rk-ksi)+cfreq(i1+1,:)*ksi

   bval=magnf(i1)*(1.0_rk-ksi)+magnf(i1+1)*ksi

END SUBROUTINE values
!--------------------



!--------------------
SUBROUTINE transf(rtor,rad,numpt,rad1,numpt1,nsort,p,powr)
!--------------------
! Transfers power array to the balance mesh
!-------
   USE numberformat
   USE gaussdata
   IMPLICIT none
   INTEGER,INTENT(IN)::numpt,numpt1,nsort
   REAL(rk),DIMENSION(numpt),INTENT(IN)::rad
   REAL(rk),DIMENSION(numpt,nsort+2),INTENT(IN)::p
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::rad1
   REAL(rk),DIMENSION(numpt1,nsort+2),INTENT(OUT)::powr
   REAL(rk):: xi,rtor,r
   INTEGER::kg,i,i1,icur

   powr=0.0_rk
   icur=2
   DO kg=1,kgmax
      r=rad1(1)+(rad1(2)-rad1(1))*ksi(kg)

      DO i=icur,numpt

         IF (r+rtor<rad(i)) EXIT

      END DO

      icur=i
      xi=(r+rtor-rad(i-1))/(rad(i)-rad(i-1))


      powr(1,:)=powr(1,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
         (1.0_rk-ksi(kg))
      powr(2,:)=powr(2,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
         ksi(kg)*(rad1(2)-rad1(1))/(rad1(3)-rad1(1))

   END DO

   DO i1=2,numpt1-2


      DO kg=1,kgmax
         r=rad1(i1)+(rad1(i1+1)-rad1(i1))*ksi(kg)

         DO i=icur,numpt

            IF (r+rtor<rad(i)) EXIT

         END DO

         icur=i
         xi=(r+rtor-rad(i-1))/(rad(i)-rad(i-1))


         powr(i1,:)=powr(i1,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
            (1.0_rk-ksi(kg))*(rad1(i1+1)-rad1(i1))/(rad1(i1+1)-rad1(i1-1))
         powr(i1+1,:)=powr(i1+1,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
            ksi(kg)*(rad1(i1+1)-rad1(i1))/(rad1(i1+2)-rad1(i1))

      END DO
   END DO

   DO kg=1,kgmax
      r=rad1(numpt1-1)+(rad1(numpt1)-rad1(numpt1-1))*ksi(kg)

      DO i=icur,numpt

         IF (r+rtor<rad(i)) EXIT

      END DO

      icur=i
      xi=(r+rtor-rad(i-1))/(rad(i)-rad(i-1))


      powr(i1,:)=powr(i1,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
         (1.0_rk-ksi(kg))*(rad1(numpt1)-rad1(numpt1-1))/(rad1(numpt1)-rad1(numpt1-2))
      powr(i1+1,:)=powr(i1+1,:)+2.0_rk*weigh(kg)*(p(i-1,:)*(1.0_rk-xi)+p(i,:)*xi)*&
         ksi(kg)

   END DO


END SUBROUTINE transf
!--------------------




!-------------------
SUBROUTINE totalpower(rad1,p,numpt1,tp)
   USE numberformat
   USE rfdata
   USE constants
   IMPLICIT NONE
   INTEGER,INTENT(IN)::numpt1
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::rad1,p
   REAL(rk),INTENT(OUT)::tp
   REAL(rk)::h
   INTEGER::i
   tp=(0.0_rk,0.0_rk)

   DO i=2,numpt1-1

      h=(rad1(i+1)-rad1(i-1))/2.0_rk

      tp=tp+4.0_rk*pi**2*rstar*(rad1(i)+rtor)*p(i)*h

   END DO

   tp=tp+4.0_rk*pi**2*rstar*(rad1(1)+rtor)*p(1)*(rad1(2)-rad1(1))/2.0_rk
   tp=tp+4.0_rk*pi**2*rstar*(rad1(numpt1)+rtor)*p(numpt1)*&
      (rad1(numpt1)-rad1(numpt1-1))/2.0_rk


END SUBROUTINE totalpower
!-------------------



!--------------------
SUBROUTINE transf1(rtor,rad,numpt,rad1,numpt1,nsort,p,powr)
!--------------------
! Transfers power array to the balance mesh
!-------
   USE numberformat
   IMPLICIT none
   INTEGER,INTENT(IN)::numpt,numpt1,nsort
   REAL(rk),DIMENSION(numpt),INTENT(IN)::rad
   REAL(rk),DIMENSION(numpt,nsort+2),INTENT(IN)::p
   REAL(rk),DIMENSION(numpt1),INTENT(IN)::rad1
   REAL(rk),DIMENSION(numpt1,nsort+2),INTENT(OUT)::powr
   REAL(rk):: ksi,rtor
   INTEGER::i,i1,icur


   powr(1,:)=p(1,:)
   powr(numpt1,:)=p(numpt,:)

   icur=2

   DO i1=2,numpt1-1

      DO i=icur,numpt

         IF (rad(i)>=rad1(i1)+rtor) EXIT

      END DO

      icur=i

      ksi=(rad1(i1)+rtor-rad(i-1))/(rad(i)-rad(i-1))


      powr(i1,:)=p(i-1,:)*(1.0_rk-ksi)+p(i,:)*ksi

   END DO

END SUBROUTINE transf1
!--------------------


!--------------------Added Tom 20190225

! This is a simple function to search for an available unit.
! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
! The UNIT value is returned by the function, and also by the optional
! argument. This allows the function to be used directly in an OPEN
! statement, and optionally save the result in a local variable.
! If no units are available, -1 is returned.

integer function newunit(unit)
   integer, intent(out), optional :: unit
! local
   integer, parameter :: LUN_MIN=10, LUN_MAX=1000
   logical :: opened
   integer :: lun
! begin
   newunit=-1
   do lun=LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
         newunit=lun
         exit
      end if
   end do
   if (present(unit)) unit=newunit
end function newunit
