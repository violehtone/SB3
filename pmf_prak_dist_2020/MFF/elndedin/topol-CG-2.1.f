c234567890

      integer i,j,k,l,tmpline,tmpint,tmpint1,tmpint2,tmpval
      integer countat,begin,endd,countchgrp
      integer atom1,atom2,atom3,atom4
      real tmpdist,tmpdistx,tmpdisty,tmpdistz
 
c DATA for the pdb file
      integer pdbnres,pdbnatom,pdbatnumb(10000),pdbresnum(1000)
      real ax(10000),ay(10000),az(10000)
      character*4 tmp
      character*3 pdbatname(10000),pdbresname(1000)
c CG FF data
      character*3 resname(20)
      integer nbeads(20),nbond(20),nangle(20),ndihe(20),nconst(20)
      character*3 tmpchar,atomname(20,7)
      character*3 beadtype(20,7)
      integer chargrp(20,7)
      real charbead(20,7)
      integer bondatom1(20,10),bondatom2(20,10),bondtype(20,10)
      real bondist(20,10),bondk(20,10)
      integer constatom1(20,10),constatom2(20,10),consttype(20,10)
      real constdist(20,10)
      integer angleatom1(20,10),angleatom2(20,10),angleatom3(20,10)
      integer angletype(20,10)
      real angle(20,10),anglek(20,10)
      integer diheatom1(20,10),diheatom2(20,10),diheatom3(20,10)
      integer diheatom4(20,10),dihetype(20,10)
      real dihe(20,10),dihek(20,10)
c DATA for topology (top)
      character*3 topatname(10000),tpresname(1000),CGname(1000)
      real atcharge(1000,7)
      integer tprestype(1000),tpnbeads(1000),beginres(1000)
      character*3 tpbeadtype(1000,7),tpatomname(1000,7)
      integer tpatnumb(1000,7),tpchrgrp(1000,7)
      real tpatcharge(1000,7)
      real*8 vector1(3),vector2(3),angleCACACA,angleCACACB,angleCBCACA
c Things for ElNeDyn
      integer elnedynkf
      real elnedyncutoff

c Paramters
      parameter     (pi=3.141592654)


C--------------------------------------------------------------
C This program is ment to read a pdb file and generate a 
C topology (protein.itp) file corresponding to the topology of 
C FFCG-2.1.
C All comments to x.periole@rug.nl; NYC-CSICUNY-03-05-2007
C--------------------------------------------------------------
c  Important Notes for Correct Functioning:
c  0- OF COURSE the use of this program is at your own risks
c  1- the pdb file given as standard imput should have been
c     generated by the program pdb2CGpdb-2.1.f
c  3- Notably the first line of the pdb file must contain the 
c     number of residues and the number of atoms in the format:
c 109   format('HEADER number of residues: ',I4," number of beads: ",I4)
C--------------------------------------------------------------

C  Define stuff for the CG model

C  Integrate the info of the CG FF
c  Note : the resname are red in cg.dat and then they are 
c  refered as the index as the order of appearance
C 

      open(unit=11,name="cg-2.1.dat",status='old')
c     read(11,*) tmpline

      do i=1,20
      read(11,222) resname(i),nbeads(i),nbond(i),nangle(i),ndihe(i)
     .,nconst(i)
c     write(6,222) resname(i),nbeads(i),nbond(i),nangle(i),ndihe(i)
c    .,nconst(i)
222   format(A3,4X,I1,3X,I1,3X,I1,3x,I1,3x,I1)
      do j=1,nbeads(i)
      read(11,*) tmpint1,beadtype(i,j),tmpint2,tmpchar,atomname(i,j)
     .,chargrp(i,j),charbead(i,j)
c     write(6,*) tmpint1,beadtype(i,j),tmpint2,tmpchar,atomname(i,j)
c    .,chargrp(i,j),charbead(i,j)
      enddo
      do j=1,nbond(i)
      read(11,*) bondatom1(i,j),bondatom2(i,j),bondtype(i,j)
     .,bondist(i,j),bondk(i,j)
c     write(6,*) bondatom1(i,j),bondatom2(i,j),bondtype(i,j)
c    ,,bondist(i,j),bondk(i,j)
      enddo
      do j=1,nconst(i)
      read(11,*) constatom1(i,j),constatom2(i,j),consttype(i,j)
     .,constdist(i,j)
c     write(6,*) constatom1(i,j),constatom2(i,j),consttype(i,j)
c    .,constdist(i,j)
      enddo
      do j=1,nangle(i)
      read(11,*) angleatom1(i,j),angleatom2(i,j),angleatom3(i,j)
     .,angletype(i,j),angle(i,j),anglek(i,j)
c     write(6,*) angleatom1(i,j),angleatom2(i,j),angleatom3(i,j)
c    .,angletype(i,j),angle(i,j),anglek(i,j)
      enddo
      do j=1,ndihe(i)
      read(11,*) diheatom1(i,j),diheatom2(i,j),diheatom3(i,j)
     .,diheatom4(i,j),dihetype(i,j),dihe(i,j),dihek(i,j)
c     write(6,*) diheatom1(i,j),diheatom2(i,j),diheatom3(i,j)
c    .,diheatom4(i,j),dihetype(i,j),dihe(i,j),dihek(i,j)
      enddo
      enddo

C  Read the pdb file -

      open(unit=13,name="pdb-CG.pdb",status='old')
      read(13,109)  pdbnres,pdbnatom
109   format('HEADER number of residues: ',I4," number of beads: ",I4)
      if (pdbnatom .le. 0 ) then 
      write(6,*) "The PDB file you give as an input is not formatted
     .correctly"
      write(6,*) "Tip: The first line should be the number of residues`
     .and then of atoms."
      stop
      endif

      do i = 1,pdbnatom
c     read (5,*) tmp,pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i),ax(i),
c    .ay(i),az(i)
      read (13,662) pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i)
     .,ax(i),ay(i),az(i)
c     write (6,662) pdbatnumb(i),pdbatname(i),pdbresname(i),pdbresnum(i),ax(i),
c    .ay(i),az(i)
      enddo

C  filling the data for each residue

      countat = 1
      countchgrp = 0
c     write(6,*) "number of resideus",pdbnres
      do i = 1,pdbnres

       beginres(i)=countat
       tpresname(i)=pdbresname(beginres(i))

        do j =1,20
        if (tpresname(i) .eq. resname(j)) then 
        tprestype(i) = j
c       write(6,*) tpresname(i),resname(j),tprestype(i),countat
        endif
        enddo
                                                                                                               
        if ((tprestype(i) .lt. 1) .and. (tprestype(i) .gt. 20)) then
        write (6,*) "Something is wrong with the attribution of the
     .residue index corresponding to its name in the pdb file!"
        stop
        endif
c                                                                                                              
c       do j =1,nbeads(tprestype(i))
         tpnbeads(i)=nbeads(tprestype(i))
         do k=1,tpnbeads(i)
          tpbeadtype(i,k)=beadtype(tprestype(i),k)
          tpatomname(i,k)=atomname(tprestype(i),k)
          tpatnumb(i,k)=countat
          tpatcharge(i,k)=charbead(tprestype(i),k)
          if ((k .gt. 1) .and. (chargrp(tprestype(i),k) .eq. 
     .       chargrp(tprestype(i),(k-1))) ) then 
           countchgrp=countchgrp
          else
           countchgrp=countchgrp+1
          endif
          tpchrgrp(i,k) = countchgrp
         countat=countat+1
c        write(6,*) i,k,tpnbeads(i),tpbeadtype(i,k),tpatomname(i,k)
c    .,tpatnumb(i,k),tpatcharge(i,k)
         enddo
c       enddo

C End of others
c     endif

      enddo

C-------------------------------------------
c  WRITING THE TOPOLOGY
C-------------------------------------------

      open(name="protein.itp",status="new",unit=12)

      write(12,666)
666   format('[ moleculetype ]')
      write(12,667)
667   format('; molname 	nrexcl')
      write(12,668)
668   format('protein      	1')
      write(12,*) ""
      write(12,669)
669   format('[ atoms ]')
      write(12,670)
670   format(';id type resnr residu atom cgnr charge')

C  Here the Atom number, type, name, charge ...

      do i = 1,pdbnres
       begin=beginres(i)
       endd=beginres(i)+nbeads(tprestype(i))-1
c      write(6,*) i,begin,endd
       do j = 1,nbeads(tprestype(i))
       tmpint2=tmpint2+chargrp(i,j)
      write(12,999) tpatnumb(i,j),tpbeadtype(i,j),i
     .,pdbresname(begin),tpatomname(i,j),tpchrgrp(i,j),tpatcharge(i,j)
       enddo
c      do j = begin,endd
c     write(12,999) pdbatnumb(j),CGname(j),pdbresnum(j),pdbresname(j)
c    .,pdbatname(j),i,atcharge(j,1)
c      enddo
      enddo

C  Here the bonds associated with the CA-CA distances

c  first the CA i-i+1

      write(12,*) ""
      write(12,998)
998   format('[ bonds ]')
      write(12,997)
997   format('; bonds between the consecutive Calphas: i-i+1')
      write(12,996)
996   format('; The distance is taken from the pdb file')

      do i = 1,(pdbnres-1)
       tmpint1=beginres(i)
       tmpint2=beginres(i+1)
       tmpdistx = ax(tmpint1)-ax(tmpint2)
       tmpdisty = ay(tmpint1)-ay(tmpint2)
       tmpdistz = az(tmpint1)-az(tmpint2)
       tmpdist = sqrt((tmpdistx**2)+(tmpdisty**2)+(tmpdistz**2))
c note that the distance is writen in nm and not Angtrom
c to follow the units used in GROMACS
      write (12,110) beginres(i),beginres(i+1),1,tmpdist/10
      enddo

c  then the CA i-i+4

      write(12,*) ""
      write(12,995)
995   format('; bonds between the consecutive Calphas: i-i+4')
      write(12,994)
994   format('; The distance is taken from the pdb file')

c Here we ask the user to define the the force constant to use
c for the CA i-i+4 interactions; can be turned off if not 
c necessary.

      write(6,*) ""
      write(6,*) "Select a value for the Ca (i-i_4) bonds."
      write(6,*) "Type 0 (zero) to exclude this type of bonds."
      write(6,*) "40000 kJ mole-1 nm-2 seems to be working fine."
      write(6,*) ""
      write(6,*) "Value:"
      read(5,*)  tmpval
c     write(6,*)  tmpval

      if (tmpval .gt. 0.0) then
      do i = 1,(pdbnres-4)
       tmpint1=beginres(i)
       tmpint2=beginres(i+4)
       tmpdistx = ax(tmpint1)-ax(tmpint2)
       tmpdisty = ay(tmpint1)-ay(tmpint2)
       tmpdistz = az(tmpint1)-az(tmpint2)
       tmpdist = sqrt((tmpdistx**2)+(tmpdisty**2)+(tmpdistz**2))
c note that the distance is writen in nm and not Angtrom
c to follow the units used in GROMACS
      write (12,111) beginres(i),beginres(i+4),1,tmpdist/10,tmpval
      enddo
      endif

c  then the CA i-i+10
 
       write(12,*) ""
       write(12,993)
 993   format('; bonds between Calphas defining the ElNeDyn')
       write(12,992)
 992   format('; bonds are taken from the pdb file between Ca')
       write(12,600)
 600   format('; at least i-i+3')

C here we read the force constant given to ElNeDyn:
c Elastic Network in Dynamics
c see: X. Periole and MA Ceruso in NSMB

       write(6,*) "Choose a force constant (kJ mol-1 nm-2) for ElNeDyn:"
       read(5,*) elnedynkf
       write(6,*) "Choose a cutoff (nm) for ElNeDyn:"
       read(5,*) elnedyncutoff
                                                                                                                
       do i = 1,pdbnres
       do j = i+3,pdbnres
c      write(12,*) "paptate"
        tmpint1=beginres(i)
        tmpint2=beginres(j)
        tmpdistx = ax(tmpint1)-ax(tmpint2)
        tmpdisty = ay(tmpint1)-ay(tmpint2)
        tmpdistz = az(tmpint1)-az(tmpint2)
        tmpdist = sqrt((tmpdistx**2)+(tmpdisty**2)+(tmpdistz**2))
        tmpdist = tmpdist/10
c note that the distance is writen in nm and not Angtrom
c to follow the units used in GROMACS
       if(tmpdist .lt. elnedyncutoff) then
       write (12,111) tmpint1,tmpint2,1,tmpdist,elnedynkf
       endif
       enddo
       enddo

c Now the bonds internal to a residue. All the bonds associated to a 
c residue are considered and then for each residue one after the other. 
c Note that the SC-SC bonds were originally constrains in the CGFF-2.0
c They are here described as bonds!

      write(12,*) ""
      write(12,991)
991   format('; bonds between Ca and SC1 for each residue')
      write(12,990)
990   format('; The distances are taken from the FF topology')

      do i = 1,(pdbnres)
       tmpint = tprestype(i)
       if (nbond(tmpint) .gt. 0) then 
        do j = 1,nbond(tmpint)
        atom1=beginres(i)+bondatom1(tmpint,j)-1
        atom2=beginres(i)+bondatom2(tmpint,j)-1
        write(12,112) atom1,atom2,bondtype(tmpint,j),
     .bondist(tmpint,j),bondk(tmpint,j)
        enddo
       endif
      enddo

      write(12,*) ""
      write(12,987)
987   format('[ constraints ]')
      write(12,989)
989   format('; constraints between atoms of the side chains')
      write(12,988)
988   format('; The values are taken from the FF topology')
                                                                                                               
      do i = 1,(pdbnres)
       tmpint = tprestype(i)
       if (nconst(tmpint) .gt. 0) then
        do j = 1,nconst(tmpint)
        atom1=beginres(i)+constatom1(tmpint,j)-1
        atom2=beginres(i)+constatom2(tmpint,j)-1
        write(12,113) atom1,atom2,consttype(tmpint,j),
     .constdist(tmpint,j)
        enddo
       endif
      enddo
                                                                                                               
      write(12,*) ""
      write(12,986)
986   format('[ angles ]')
      write(12,985)
985   format('; angles between atoms of the side chains')
      write(12,984)
984   format('; The angles are taken from the FF topology')
                                                                                                               
      do i = 1,(pdbnres)
       tmpint = tprestype(i)
       if (nangle(tmpint) .gt. 0) then
        do j = 1,nangle(tmpint)
        atom1=beginres(i)+angleatom1(tmpint,j)-1
        atom2=beginres(i)+angleatom2(tmpint,j)-1
        atom3=beginres(i)+angleatom3(tmpint,j)-1
        write(12,114) atom1,atom2,atom3,angletype(tmpint,j),
     .angle(tmpint,j),anglek(tmpint,j)
        enddo
       endif
      enddo

      write(12,*) ""
      write(12,983)
983   format('; Here the angles between CAs i-1/i/i+1')
      write(12,982)
982   format('; The values are taken from the struture')
                                                                                                               
      do i = 2,(pdbnres-1)

c  Note that the angles are centered on the atom1, here i
c  And written i-1, i, i+1

      atom1=beginres(i)
      atom2=beginres(i-1)
      atom3=beginres(i+1)

c--- there the vector Ca i-1 -> i
      vector1(1) = ax(atom2)-ax(atom1)
      vector1(2) = ay(atom2)-ay(atom1)
      vector1(3) = az(atom2)-az(atom1)
       
c--- there the vector Ca i -> i+1
      vector2(1) = ax(atom3)-ax(atom1)
      vector2(2) = ay(atom3)-ay(atom1)
      vector2(3) = az(atom3)-az(atom1)

c---- now the angle calculation

      angleCACACA=angvec(vector1,vector2)
      angleCACACA=angleCACACA*180.0/pi

c---- now writing the angles

      write(12,116) atom2,atom1,atom3,2,angleCACACA,40.0
      write(12,116) atom2,atom1,atom3,2,120.0,40.0

      enddo
                                                                                                               
C Not these anymore !
c      write(12,*) ""
c      write(12,977)
c977   format('; Here the angles between CA(i-1)-CA(i)-CB(i)')
c      write(12,976)
c976   format('; The values are taken from the struture')
c      write(12,975)
c975   format('; Note that CB is the second bead of the SC')
c      write(12,974)
c974   format('; NOT necesseraly the CB in the structure')
c
c      do i = 2,(pdbnres)
c
c      tmpint = tprestype(i)
c
cc--- this is to exclude the GLY and ALA which have no CB
c      if ((resname(tmpint) .ne. "GLY") .and. 
c     .(resname(tmpint) .ne. "ALA")) then 
c
cc  Note that the angles are centered on the atom 1, here i
cc  And written i-1, i, (i)+1
c
c      atom1=beginres(i)
c      atom2=beginres(i-1)
c      atom3=beginres(i)+1
c
cc--- there the vector Ca i-1 -> i
c      vector1(1) = ax(atom2)-ax(atom1)
c      vector1(2) = ay(atom2)-ay(atom1)
c      vector1(3) = az(atom2)-az(atom1)
c
cc--- there the vector Ca i -> Cb i
c      vector2(1) = ax(atom3)-ax(atom1)
c      vector2(2) = ay(atom3)-ay(atom1)
c      vector2(3) = az(atom3)-az(atom1)
c
cc---- now the angle calculation
c
c      angleCACACB=angvec(vector1,vector2)
c      angleCACACB=angleCACACB*180.0/pi
c
cc---- now writing the angles
c
c      write(12,117) atom2,atom1,atom3,2,angleCACACB,10.0
c
c      endif
c
c      enddo
c
c      write(12,*) ""
c      write(12,973)
c973   format('; Here the angles between CB(i)-CA(i)-CA(i+1)')
c      write(12,972)
c972   format('; The values are taken from the struture')
c      write(12,971)
c971   format('; Note that CB is the second bead of the SC')
c      write(12,970)
c970   format('; NOT necesseraly the CB in the structure')
c
c      do i = 1,(pdbnres-1)
c
c      tmpint = tprestype(i)
c
cc--- this is to exclude the GLY and ALA which have no CB
c      if ((resname(tmpint) .ne. "GLY") .and.
c     .(resname(tmpint) .ne. "ALA")) then
c 
cc  Note that the angles are centered on the atom1, here i
cc  And written (i)+1, i, i+1
c
c      atom1=beginres(i)
c      atom2=beginres(i)+1
c      atom3=beginres(i+1)
c
cc--- there the vector Cb i -> Ca i
c      vector1(1) = ax(atom2)-ax(atom1)
c      vector1(2) = ay(atom2)-ay(atom1)
c      vector1(3) = az(atom2)-az(atom1)
c
cc--- there the vector Ca i -> i+1
c      vector2(1) = ax(atom3)-ax(atom1)
c      vector2(2) = ay(atom3)-ay(atom1)
c      vector2(3) = az(atom3)-az(atom1)
c
cc---- now the angle calculation
c
c      angleCBCACA=angvec(vector1,vector2)
c      angleCBCACA=angleCBCACA*180.0/pi
c
cc---- now writing the angles
c
c      write(12,117) atom2,atom1,atom3,2,angleCBCACA,10.0
c
c      endif

c      enddo

      write(12,*) ""
      write(12,981)
981   format('[ dihedrals ]')
      write(12,980)
980   format('; Dihedrals angles used in Side chains')
      write(12,979)
979   format('; The angles are taken from the FF topology')
                                                                                                               
      do i = 1,(pdbnres)
       tmpint = tprestype(i)
       if (ndihe(tmpint) .gt. 0) then
        do j = 1,ndihe(tmpint)
        atom1=beginres(i)+diheatom1(tmpint,j)-1
        atom2=beginres(i)+diheatom2(tmpint,j)-1
        atom3=beginres(i)+diheatom3(tmpint,j)-1
        atom4=beginres(i)+diheatom4(tmpint,j)-1
        write(12,115) atom1,atom2,atom3,atom4,dihetype(tmpint,j),
     .dihe(tmpint,j),dihek(tmpint,j)
        enddo
       endif
      enddo


      write(12,*) ""
      write(12,108)
108   format('#ifdef POSRES')
      write(12,107)
107   format('#include "posre.itp"')
      write(12,106)
106   format('#endif')

110   format(I4,4X,I4,4X,I4,4X,F8.3,4X,'150000')
111   format(I4,4X,I4,4X,I4,4X,F8.3,4X,I8)
112   format(I4,4X,I4,4X,I4,4X,F8.3,4X,F8.1)
113   format(I4,4X,I4,4X,I4,4X,F8.3)
114   format(I4,4X,I4,4X,I4,4X,I4,4X,F8.1,4X,F8.1)
115   format(I4,4X,I4,4X,I4,4X,I4,4X,I4,4X,F8.1,4X,F8.1)
116   format(I4,4X,I4,4X,I4,4X,I4,4X,F8.1,4X,F8.1)
117   format(I4,4X,I4,4X,I4,4X,I4,4X,F8.1,4X,F8.1)
999   format(I4,4X,A3,4X,I4,5X,A3,4X,A3,I4,6X,F5.2)
662   format('ATOM',I7,1X,A4,1X,A3,2X,I4,4X,3F8.3)

      end

C-----Here some usefull functions for the angles

      subroutine prodvec(v1,v2,v3)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      implicit none
      real*8   v1(3), v2(3),v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      end

      function prodsca(v1,v2)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      implicit none
      real*8 prodsca
      real*8 v1(3),v2(3)

      prodsca=0.0D+00
      prodsca=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      end

      function vecnorm(v)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      implicit none
      integer i
      real*8 vecnorm,v(3)

      vecnorm=0.0D+00

      do i=1,3
        vecnorm= vecnorm + v(i)*v(i)
      enddo

      vecnorm=dsqrt(vecnorm)
      end

      function angvec(v1,v2)
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      implicit none
      integer i
      real*8 angvec
      real*8 v1(3),v2(3),vecnorm,prodsca

      if (vecnorm(v1) .gt. 0.0 .and. vecnorm(v2) .gt. 0.0) then
        angvec=prodsca(v1,v2)/(vecnorm(v1)*vecnorm(v2))
        if (angvec.ge.1.0) then
          angvec=acos(1.0)
        else if (angvec.le.-1.0) then
          angvec=acos(-1.0)
        else
          angvec=acos(angvec)
        end if
      else
        write(0,*) 'Vector with zero norm found ???'
        stop
      end if

      end