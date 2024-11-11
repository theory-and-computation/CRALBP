      program wham1sp   
c
c This program performs one-dimensional Weighted Histogram Amplitude Method
c to unbias the windows of an umbrella sampling calculation with harmonic functions
C Version WHAM1SP write out a file of the F(i)
C Use symmetrized histograms if periodicity different from zero                  
c
C Compile with f77 -O2 -o wham1sp wham1sp.f -lm 
C Syntax ./wham1sp <toto.list >toto.w1sp   
C
C Changed "open(unit=XXX,file=fnam,status='unknown')"
C      to "open(unit=XXX,name=fnam,status='unknown')"
C   to be compatible with gfortran.  -Albert Lau
C Compile with gfortran wham1sp.f -o wham1sp
C
C Example of input file: toto.list 
C
C TOTO2P.PMF   Output PMF file name 
C TOTO2P.RHO   Output unbiased distribution file name 
C TOTO2P.BIA   Output symmetrized biased distribution file name
C TOTO.FFF     File containing the Fis
C   -190.0 190.0  10.0  360.0  xmin, xmax, delta, Period (P1)
C                              Set the period to 0 for non periodic WHAM  
C   15 1000 100     Number of windows, max number of iterations, how often to 
C                   write the Fi's
C  0.001            Tolerance  
C PMFBR9_1.TIM      First time serie 
C -180  1.218468E-02   x1, k1 
C PMFBR9_2.TIM
C -150  1.218468E-02 
C .....................................
C PMFBR9_15.TIM     Last time series 
C -180  1.218468E-02  
C
C Example of time series file
C# 1500 STEPS BACT, DIHE1=187.5  FORCE1=120.0   Comment line(s)
C  0.0010  183.493554    <time> x1           (time is optional) 
C  0.0020  183.335161   
C-----------------------------------
C  1.4920  186.834509  
C
C The program bins the time series into 1-D histograms for the variable
C x1, extending from xmin to xmax 
C Points outside this range are discarded (P=0) or symmetrized (P != 0).
C Bin widths is delta -> number of bins = (xmax-xmin)/delta
C This number of bins is supposed to be less than the adjustable parameter MAX2 
C MAX1 is the maximum number of allowed time-series.
C Each time serie  results from an umbrella sampling centered at (x1) 
C with force constant  k1 
C Biasing potential k1(x-x1)**2 ... (Careful : No factor 0.5) ... 
C If several time series are entered with same force constants and centered 
C on the same value of the reaction coordinate modulo P1, the 
C corresponding histogramms will be combined automatically.
C The time series may be given in any order in the input file 
C The program exits when whichever comes first : tolerance or 
C max. number of iterations is reached 
C If the number of counts in one bin of the biased histogram is zero, then the
C PMF is undefined for that bin and set to an arbitrarily high value
C namely : 9999.999
C
      IMPLICIT NONE
      INTEGER MAX1, MAX2 
C IF YOU CHANGE THESE PARAMETERS DON'T FORGET TO CHANGE THEM LINE 343 ALSO ! 
      parameter(MAX1=1000,MAX2=4001)
      CHARACTER*80 fnam
      CHARACTER*80 tseries(MAX1)
      INTEGER H(MAX1, MAX2) 
      INTEGER map(MAX2)
      CHARACTER*20 hisfor, fifor  
      logical conv
      INTEGER nb_data(MAX1)
      INTEGER ihist, Nw, nbin1
      INTEGER   j, jmin
      INTEGER i, icycle, tbin1
      INTEGER NIter , fifreq 
      INTEGER ix1 

      DOUBLE PRECISION    rho(MAX2), pmf(MAX2)
      DOUBLE PRECISION    Top(MAX2)
      DOUBLE PRECISION    V1(MAX2,MAX1)
      DOUBLE PRECISION    x1(MAX2)
      DOUBLE PRECISION    F(MAX1), F2(MAX1)
      DOUBLE PRECISION    V(MAX1)
      DOUBLE PRECISION    k1(MAX1), cx1(MAX1)
      DOUBLE PRECISION    xmin, xmax, delta
      DOUBLE PRECISION    x1m
      DOUBLE PRECISION    rhomax, pmfmin
      DOUBLE PRECISION    bot, sum, Diff, Diff0, Difft
      DOUBLE PRECISION    tol
      DOUBLE PRECISION    kboltz, kbt, temp
      DOUBLE PRECISION    kbtm1
      DOUBLE PRECISION    P1, halfP1, xmcx1  
C  Default temperature 300 Kelvin 
      parameter(kboltz=1.9872D-3,temp=300.0D0)
C     Common block
      common /histo/ H

c Boltzmann factor
      kbt=kboltz*temp   
      kbtm1=1.0D0/kbt
      tol = 0.00002
c
      write(*,*) '--------------------------------------------------'
      write(*,*) 'ONE-DIMENSIONAL UMBRELLA SAMPLING ANALYSIS PROGRAM'
      write(*,*) 'Authors:  Serge Crouzy and Benoit Roux from 2-D '
      write(*,*) 'version by Benoit Roux (1994).'
      write(*,*) 'Periodic version: Serge Crouzy (1999-2001).'
      write(*,*) '--------------------------------------------------'
      write(*,*)
c

      read(*,'(a)') fnam
      call lower(fnam)
      open(unit=20,file=fnam,status='unknown')
      write(*,*) 'Free energy profile output files: ',fnam(1:25)

      read(*,'(a)') fnam
      call lower(fnam)
      open(unit=21,file=fnam,status='unknown')
      write(*,*) 'Average density output files: ',fnam(1:25)

      read(*,'(a)') fnam
      call lower(fnam)
      open(unit=23,file=fnam,status='unknown')
      write(*,*) 'Biased distribution output file: ',fnam(1:25)

      READ(*,'(a)') fnam
      call lower(fnam)
      OPEN(unit=27,file=fnam,status='unknown')
      WRITE(6,*) 'F(i) output file: ',fnam(1:25)

      read(*,*) xmin,xmax,delta,P1 
c
      write(*,*) 'The absolute minimum of the PMF will be set to zero'

      write(*,*)
      write(*,100) 'first position in the histogram ', xmin
      write(*,100) 'last position in the histogram ', xmax
      write(*,100) 'Period ',P1
      write(*,100) 'width of the bins       ', delta
      nbin1=int((xmax-xmin+0.5*delta)/delta)
      write(*,101) 'number of bins                    ', nbin1
      write(*,*)

  100 format(1x,a,f8.3)
  101 format(1x,a,i8)

      write(*,*)
      if (nbin1.gt.MAX2) then
      write(*,*) 'Number of bins too large'
      stop
      ENDIF
      WRITE(6,*) '--------------------------------------------------'
C
      read(*,*) Nw, NIter, fifreq 
      write(*,101) 'Total number of windows           ', Nw
      if (Nw.gt.MAX1) then
      write(*,*) 'number of histograms too large'
      stop
      ENDIF
      WRITE(6,'(/,1x,a,i8,a)') 'Save the F(i) every ',fifreq,' steps. '
      WRITE(fifor,136) Nw
 136  FORMAT('(I5,',I3,'F8.3)')
      write(*,101) 'Number of iterations requested    ', NIter
C     read tolerance
      read(*,*) tol
      write(6,'(1x,a,F8.5)') 'Tolerance for PMF calculation ',tol
      WRITE(6,*) '--------------------------------------------------'

C     Initialize the biased distribution
      do 203 ix1=1,MAX2
          map(ix1)=0
 203  CONTINUE
      write(hisfor,135) nbin1
 135  format('(',I9,'I9)')

c Read the list of time series
      write(*,*) 'Read the list of time series'
c open the files of the time series time-series
      do 764 ihist=1,Nw
      F2(ihist) = 0.0D0
      read(*,'(a)') fnam
      call lower(fnam)
      tseries(ihist)=fnam
      read(*,*)  cx1(ihist),k1(ihist)
      write(*,'(1x,a)') tseries(ihist)(1:30)
      write(*,'(1x,a,i5,a,4f10.3)') 'Parameters of window potential',
     ,     ihist,':',k1(ihist),cx1(ihist)
 
c Make the periodic histograms
      open(unit=10,file=tseries(ihist),status='old')
      call mkH(10,ihist,nb_data, 
     ,     delta,xmin,xmax,P1,nbin1)
      close(unit=11)
  764 CONTINUE
      close(unit=31)
      WRITE(6,*) '--------------------------------------------------'
C
C     Calculate and Print the symmetrized biased distribution
      DO i=1,Nw
        DO j=1,nbin1
            map(j)=map(j)+H(i,j)
        ENDDO
      ENDDO

      do 207 ix1=1,nbin1
          write(23,hisfor) map(ix1)
 207  CONTINUE
      close(unit=23)
C
c Iterate the WHAM equation to unbias and recombine the histograms
C Some calculations valid for all cycles
      do 211 j=1,nbin1

      x1(j)=xmin+j*delta-0.5*delta
c     write(*,*) 'x1 ',x1(j)

      Top(j) = 0.0D0

      do 212 ihist=1,Nw
        Top(j) = Top(j) + DBLE(H(ihist,j))
  212 CONTINUE
  211 CONTINUE

C Here we must make an assumption : we assume that all histograms
C have empty bins outside [center-P/2,center+P/2]
C
      halfP1 = 0.5 * P1
      do 311 ihist=1,Nw
      do 312 j=1,nbin1
      xmcx1 = DABS(x1(j)-cx1(ihist))
      if (xmcx1.gt.halfP1) xmcx1 = xmcx1 - P1
      V1(j,ihist) = k1(ihist)*xmcx1*xmcx1
  312 CONTINUE
  311 CONTINUE

      write(*,*) 'Begin to iterate on the F()'
      icycle = 0
      Difft = 1.0D10

 1000 icycle=icycle+1

      do 202 ihist=1,Nw
        F(ihist) = F2(ihist)
        F2(ihist) = 0.0D0
  202 CONTINUE 

      do 21 j=1,nbin1

      Bot = 0.0D0

      do 22 ihist=1,Nw
        V(ihist) = V1(j,ihist)
        
        Bot = Bot + 
     ,    DBLE(nb_data(ihist)) * DEXP((F(ihist)-V(ihist))*kbtm1)
   22 CONTINUE 

      rho(j) = Top(j)/Bot
C      write(*,*) 'j, rho(j)',j,rho(j)
C      write(*,*) 'Top(j), Bot',Top(j), Bot   
c
c Now sum over rbins to get updated free energies F(ihist)
c
      do 23 ihist=1,Nw
        F2(ihist) = F2(ihist) + rho(j)*DEXP(-V(ihist)*kbtm1)
C        write(6,*) 'ihist F2(ihist) ',ihist,F2(ihist) 
   23 CONTINUE 
   21 CONTINUE
C      write(*,*) Top(1), rho(1) 
C      write(*,*) Top(nbin1), rho(nbin1) 
c 
c test for convergence 
c
      conv = .true.
      Diff0=0.0D0
      do 24 ihist=1,Nw
        F2(ihist) = -kbt*DLOG( F2(ihist) )
        Diff = Abs( F2(ihist) - F(ihist) )
        if (Diff.gt.tol) Conv = .False.
        if (Diff .gt. Diff0) Diff0=Diff
C      write(*,*) 'Cycle=',icycle,' Window=',ihist,
C     ,            ' F=',(F(ihist)-F(1))
   24 CONTINUE
      Write(6,*) 'Round = ',icycle,'  Diff = ',Diff0
      if (Diff0 .lt. Difft) then
        Difft = Diff0
      else
        Difft = Diff0
        write(6,*) 'Difference has increased since previous step
     + - Go on anyway !'
      end if
      if (fifreq .ne.0 .and. (mod(icycle,fifreq) .eq. 0 .or. 
     + icycle .eq. NIter .or. Conv)) THEN
        WRITE(27,fifor) icycle,(F2(ihist), ihist=1,Nw)
      ENDIF
      if ((icycle.lt.NIter).and.(.Not.Conv)) Goto 1000
c
      write(6,*) 'Number of cycles:',icycle
      do 25 ihist=1,Nw
      write(6,*) 'The converged F in window ',
     ,            ihist,' = ',(F(ihist)-F(1))
   25 CONTINUE
      write(6,*)

c -------------------------------------------------------------------------
c find the maximum of the rho
      rhomax = rho(1)
      jmin = 1
      do 30 j=1,nbin1
      if (rhomax.lt.rho(j)) then
      rhomax = rho(j)
      jmin = j
      ENDIF
   30 CONTINUE
      x1m=xmin+jmin*delta-delta/2
      write(*,*) 'Maximum density at: x1 = ',x1m
      write(*,*)

      write(*,*) '----------------------------------'
      write(*,*) '# One-dimensional PMF(X1)'
      write(*,*) 'PMF set to zero for maximum density...' 

c Make the PMF from the rho

      do 26 j=1,nbin1
      if (rho(j).ne.0.0)then
      pmf(j)=-kbt*DLOG(rho(j)/rhomax)
      else
C An arbitrarily high value 
      pmf(j)=9999.999
      ENDIF
   26 CONTINUE

      do 40 j=1,nbin1
      write(20,'(6f12.3)') x1(j),pmf(j)
      write(21,'(2f12.3,f12.8)') x1(j),rho(j)/rhomax
   40 CONTINUE

      end

      subroutine mkH(iunit,ihist,nb_data, 
     ,           delta,xmin,xmax,P1,nbin1)

c  This program makes histograms from files coming out of CORREL
c  written with the write DUMB TIME options
c
      INTEGER iunit,MAX1,MAX2,ihist,nread
      parameter(MAX1=1000,MAX2=4001)
      INTEGER  H(MAX1, MAX2)
      INTEGER nb_data(MAX1) 
      DOUBLE PRECISION delta,xmin,xmax
      DOUBLE PRECISION P1 
      DOUBLE PRECISION junk 
      INTEGER nbin1,iread 
      DOUBLE PRECISION    kbtm1

c Local variables
      DOUBLE PRECISION time,x1
      INTEGER sym(MAX1)
      INTEGER ix1,nb1 
      CHARACTER*85 line
      DOUBLE PRECISION sx1,sxx1
      DOUBLE PRECISION minx1, maxx1

C     Common block
      common /histo/ H

      nb1 = int((P1+0.5*delta)/delta)
      minx1=1.0D10
      maxx1=-1.0D10
c Read title if there is one
  999 read(iunit,'(a)') line
      if ((line(1:1).eq.'#').or.(line(1:1).eq.'*'))then
      write(*,*) line(1:79)
      goto 999
      else
      backspace(unit=iunit)
      ENDIF

      sx1=0.0
      sxx1=0.0
      nb_data(ihist)=0
      sym(ihist)=0 
      do 111 j=1,nbin1
            H(ihist,j)=0
  111 CONTINUE

c Check how many columns of numbers are in the file
      read(iunit,'(a)') line
      backspace(unit=iunit)
      ncol=0
      do 2 i=1,84
      if ((line(i:i).ne.' ').and.(line(i+1:i+1).eq.' '))then
      ncol=ncol+1
      ENDIF
    2 CONTINUE
      write(*,*) 'there are ',ncol,' columns of data'
      time=0.0
      x1=0.0

C      write(6,*) 'xmin ',xmin
C      write(6,*) 'MAX2 ',MAX2
C      write(6,*) 'delta ',delta
      iread=1
 1000 if (ncol.ge.3)then
C Column 3 is not used anymore 
      read(iunit,*,end=2000) time,x1,junk      
      elseif (ncol.eq.2)then
      read(iunit,*,end=2000) time,x1  
      elseif (ncol.eq.1)then
      time=9999.99 
      read(iunit,*,end=2000) x1   
      ENDIF
C      write(6,*) 'time, x1 ',time,x1
      if (x1.lt.minx1) minx1=x1
      if (x1.gt.maxx1) maxx1=x1
C   Statistics on histogram ... 
      if (P1.EQ.0.0) THEN
        if (x1.GE.xmin .AND. x1.LE.xmax) THEN
         sx1=sx1+x1
         sxx1=sxx1+x1*x1
         ix1=int((x1-xmin)/delta)+1
         IF(ix1.eq.nbin1+1) ix1 = nbin1
         nb_data(ihist)=nb_data(ihist)+1
         H(ihist,ix1)=H(ihist,ix1)+1
        ENDIF
      ELSE
        sx1=sx1+x1
        sxx1=sxx1+x1*x1
        nb_data(ihist)=nb_data(ihist)+1
        if (x1.lt.xmin) x1=x1+P1
        if (x1.gt.xmax) x1=x1-P1
        ix1=int((x1-xmin)/delta)+1
        H(ihist,ix1)=H(ihist,ix1)+1
C      write(6,*) 'ix1 ',ix1
      ENDIF
C 
C We now make the histograms completely symmetrical above P/2 and
C below -P/2.  This section is relevant if P1 < (xmax-xmin) <= 2*P1.
C
      if (P1.NE.0.0) THEN 
         if ((ix1+nb1).le.nbin1) then
            H(ihist,ix1+nb1)=H(ihist,ix1+nb1)+1
            sym(ihist)=sym(ihist)+1
         endif 
         if ((ix1-nb1).ge.1) then
            H(ihist,ix1-nb1)=H(ihist,ix1-nb1)+1
            sym(ihist)=sym(ihist)+1
         endif 
      ENDIF
 
      iread=iread+1
      goto 1000

 2000 nread=iread-1
      sx1=sx1/nb_data(ihist)
      sxx1=sxx1/nb_data(ihist)
      sxx1=sqrt(sxx1-sx1*sx1)
      write(6,*) 'Statistics for time series #',ihist  
      write(6,'(a,F10.3,a,F10.3)') ' Minx1 ',minx1,' Maxx1 ',maxx1
      write(6,*) 'Statistics for histogram #',ihist   
      write(6,'(i9,a)')  nb_data(ihist),' points in the histogram'
      write(6,102) 'Average x ', sx1,' rms     x ', sxx1
      if (P1.NE.0.0) THEN
        nb_data(ihist)=nb_data(ihist)+sym(ihist) 
        write(6,'(i7,a)')  nb_data(ihist),' points in the periodic 
     ,histogram'
  102 format(2(1x,a,f8.3))
      ENDIF 
      write(*,*)

      RETURN
      end

c------------------------------------------------------------------------------
      subroutine lower(c1)
c     Lower case substitutions.
c     The fortran function ICHAR(a) RETURNs the ascii number corresponding to 
c     the CHARACTER (a).  The function CHAR(i) RETURNs the CHARACTER 
c     corresponding to the asci number i.
c     For alphabetic CHARACTERs:
c     A = 65, B = 66, ..., Z = 90
c     a = 97, b = 98, ..., z = 122
      CHARACTER*(*) c1
      INTEGER lenc1,ic1

      lenc1=len(c1)
      do 431 i=1,lenc1
      ic1=ichar(c1(i:i))
      if ((ic1.ge.65).and.(ic1.le.90))then
      c1(i:i)=char(ic1+32)
      ENDIF
  431 CONTINUE
      RETURN
      end

