	   program  grain_construction

	   implicit real*8 (a-h,o-z)

       parameter (ngus=10000,nly=1,npos=7500)

	   parameter(nsc=npos**2)
!c	ngus:: initial guess of number of grain in each layer
!c	nly: number of scanned layers
!c	nstpsz:: is resolution of the system and also the mesh size of 
!c	each grain
!c	x0,z0,xfin,and zfin are the initial and final grid points 
!c	npos:: dimension of the pos variable
      
	   dimension grain(ngus,nly,50),ginfo(50),x1d(npos),y1d(npos), &
	   z1d(npos),grain_diam(ngus,nly),nogly(nly),pos(nsc,nly,4),&
	   eul(3),ngb(ngus),ngbgr(ngus,1000)
       
       character*256 OUTDIR


!c	graion(i,j,k):: i is number, j is layer and k is the information
!c	                 k=43 is the grain diameter
!c-----nogly(i):: number of counted grain in each layer.
!c	pos(i,j,k):: i is the counter; j is the layer; and k =1:4 is x,y,z and grain number 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       Write (*,*) "please enter step size"
       read (*,*) nstpsz
       write (*,*) "this is what I read", nstpsz
      
       Write (*,*) "please enter x0"
       read (*,*) x0
       write (*,*) "this is what I read", x0
       Write (*,*) "please enter y0"
       read (*,*) y0
       write (*,*) "this is what I read", y0
       Write (*,*) "please enter z0"
       read (*,*) z0
       write (*,*) "this is what I read", z0
       Write (*,*) "please enter xfinal"
       read (*,*) xfin
       write (*,*) "this is what I read", xfin
       Write (*,*) "please enter yfinal"
       read (*,*) yfin
       write (*,*) "this is what I read", yfin
       Write (*,*) "please enter zfinal"
       read (*,*) zfin
       write (*,*) "this is what I read", zfin
!!!!!!!!!!!!!!!!!!!!!!
       Write (*,*) "Voronoi (0) or Weighted Voronoi (1)?"
       read (*,*) nvorwvor
       if (nvorwvor .eq. 0) then
            write (*,*) "Just Voronoi"
       else 
           write (*,*) "Weighted Voronoi"
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       Write (*,*) "please enter the name of the file"
       read (*,*) OUTDIR
       OUTDIR=trim(OUTDIR)
       write (*,*) "this is what I read  ", trim(OUTDIR)
      
      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
       
!	   open(UNIT=1,FILE='esrf5n.dat')
	   open(UNIT=1,FILE= trim(OUTDIR))
       open(UNIT=13,FILE='boundaries.dat')
	   open(UNIT=40, FILE='neighbours.dat')
	


100	   format(42(f16.5))
107	   format(100(2x,f14.5))
108	   format(200(2x,I14))



	   do i=1,nly
	     nogly(i)=0
	
	   end do	

	   pi=4.d0*atan(1.d0)


      
      
      
!c-----reading info from external file
	   do j=1,nly
	
	     test=1
	     i=0
	     do while (test .eq. 1)
	       i=i+1
	       read(j,100,iostat=ktst) (ginfo(jj), jj=1,42)
	  
	       if (ktst .ne. 0) then
	        test=0
	        go to 200
	       else
	        nogly(j)=nogly(j)+1
	       do k=1,42
			if (k .ge. 4 .and. k .le.6) then
	          grain(i,j,k)=ginfo(k)*1000.d0 !to convert milimeters to to micron
	        else
		      grain(i,j,k)=ginfo(k)
	        end if
	       end do
	      
	       call rodrig2euler(grain(i,j,7:9),eul)  !to convert rodrigues to Euler angles
           !!!!! Rodrigues are from local 2 global, here I convert them from global 2 local
           !!!! This is to make it consistenet with CPFE modeling and Bunge notation
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           p1=-eul(3)+360
           p=-eul(2)+360
           p2=-eul(1)+360
           
           eul(1)=p1
           eul(2)=p
           eul(3)=p2
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	       do k =1,3
	        grain(i,j,k+6)=eul(k) !euler angles are form sample2crystal
	       end do
	
	      end if
	    
!c	    dg= (4.d0*ginfo(3)/pi)**(1./2.) ! this is for layer wise-considering area!!!
!c	    grain(i,j,43)=dg
!          volum=0.4446  !total volume mm3
!          vr=grain(i,j,3)/33243.8
!          vi=vr*volum*1e9
!          ginfo(3)=vi

!	      dg= (6.d0*ginfo(3)/pi)**(1./3.) ! this is for 3d constructed layers-considering volume!!!
!	      grain(i,j,43)=dg
!	      grain(i,j,43)=0.d0



200	      if (ktst .ne. 0) then
	       write(*,*) '#grains= ', i-1
	       write(*,*) 'Layer=', j
	      end if
	    end do

       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This is to normalize the measured volume        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
       ngrni=i-1
       j=1
       boxvol=0.d0
       nboxngrn=0
       do i1=1, ngrni
           !!! this is to check if the grain falls into the selected volume
           if (grain(i1,j,4) .ge. x0 .and. grain(i1,j,4) .le. xfin) then
               if (grain(i1,j,5) .ge. y0 .and. grain(i1,j,5) .le. yfin) then
                   if (grain(i1,j,6) .ge. z0 .and. grain(i1,j,6) .le. zfin) then
                       nboxngrn=nboxngrn+1
                       boxvol=boxvol+grain(i1,j,3)
                                             
                   end if 
               end if
           end if 
       end do
       
       write (*,*) "number of grains counted in the box are:  ", nboxngrn
       !!! the physical diemnsion of the box 
       fdob=(xfin-x0)*(yfin-y0)*(zfin-z0)
       j=1
       do i1=1,ngrni
           vr=fdob*(grain(i1,j,3)/boxvol)
           dg=(6.d0*vr/pi)**(1./3.) 
           if (nvorwvor .eq. 0) then
               grain(i1,j,43)=0.d0
           else
               grain(i1,j,43)=dg
           end if               
       end do
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!c--------end reading and assigning data
!c---- generating a 3d matrix of points in cartesian space
	   i1d=0
	   do i=x0,xfin,nstpsz
	  
	    i1d=i1d+1
	    x1d(i1d)=i
	  	
	   end do



	   i2d=0
	   do i=y0,yfin,nstpsz
	  
	    i2d=i2d+1
	    y1d(i2d)=i
	  	
	   end do

	
	   write(*,*) 'i1d',i1d


!c----------for the 3rd direction

	   i3d=0
	   do i=z0,zfin,nstpsz
	  
	    i3d=i3d+1
	    z1d(i3d)=i
	  	
	   end do


!c--------------
	
	
	   do i=1,nly
	     ixd=0
	     do j=x0,xfin, nstpsz
	    		
	       do k=y0,yfin,nstpsz

	        do iz=z0,zfin,nstpsz
	      
		     ixd=ixd+1

		     pos(ixd,i,1)=j !x 
	         pos(ixd,i,2)=k ! y
	         pos(ixd,i,3)=iz !z
	         pos(ixd,i,4)=1   !grain #
!c	         write(*,*) j,k

	        end do
	       end do

	     end do
	   end do

	   write(*,*) 'ixd',ixd
!c
!c	write(*,*) 'izd',izd

	   do i=1,nly
	    write(*,*) "LAYER #", i
	    call grainmaker(grain,nogly,pos,i,ngus,nly,nsc,i1d,i2d,i3d)
	    end do

!c-----------------------This is to generate neighbours	
!c----------------------assuming there is jus one layer (all of the info in one layer including Z info)
	   write (*,*) 'Neighbours creation is initiated'

	   do i=1,ngus
	    ngb(i)=0 !number of neighbours
	    do j=1,1000
	     ngbgr(i,j)=0
	    end do
	   end do
	
	
	   ixd1=0
	   npx=floor((xfin-x0)/nstpsz)+1
	   npy=floor((yfin-y0)/nstpsz)+1
	   npz=floor((zfin-z0)/nstpsz)+1
	   do j=x0,xfin, nstpsz
	    		
	    do k=y0,yfin,nstpsz

	     do iz=z0,zfin,nstpsz
	      ixd1=ixd1+1
!c--------------negatiove x
	      xtry=pos(ixd1,1,1)
	      if (xtry .gt. x0) then
	        xtry=xtry-nstpsz
	        idtry=ixd1-npz*npy
	        if(pos(idtry,1,1) .ne. xtry) then
	          write(*,*) 'problem in negative x'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	         ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	         ! this is for common boundary
	         
	        end if
          end if
!c----------------positive x
	      xtry=pos(ixd1,1,1)
	      if (xtry .lt. xfin) then
	        xtry=xtry+nstpsz
	        idtry=ixd1+npz*npy
	        if(pos(idtry,1,1) .ne. xtry) then
	          write(*,*) 'problem in positive x'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	         ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	        end if
          end if	        
!c----------------------negative y
	      ytry=pos(ixd1,1,2)
	      if (ytry .gt. y0) then
	        ytry=ytry-nstpsz
	        idtry=ixd1-npz
	        if(pos(idtry,1,2) .ne. ytry) then
	          write(*,*) 'problem in negative y'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	         ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	        end if
          end if	  
!c----------------------positive y
	      ytry=pos(ixd1,1,2)
	      if (ytry .lt. yfin) then
	        ytry=ytry+nstpsz
	        idtry=ixd1+npz
	        if(pos(idtry,1,2) .ne. ytry) then
	          write(*,*) 'problem in positive y'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	         ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	        end if
          end if
!c----------------------negative z
	      ztry=pos(ixd1,1,3)
	      if (ztry .gt. z0) then
	        ztry=ztry-nstpsz
	        idtry=ixd1-1
	        if(pos(idtry,1,3) .ne. ztry) then
	          write(*,*) 'problem in negative z'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	         ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	        end if
          end if
!c----------------------positive z
	      ztry=pos(ixd1,1,3)
	      if (ztry .lt. zfin) then
	        ztry=ztry+nstpsz
	        idtry=ixd1+1
	        if(pos(idtry,1,3) .ne. ztry) then
	          write(*,*) 'problem in positive z'
	          stop
	        end if
	        if (pos(ixd1,1,4) .ne. pos(idtry,1,4)) then !not same grain
	       	 ncrg=pos(ixd1,1,4)
	         ngcr=pos(idtry,1,4)
	         write(13,107) float(ncrg),float(ngcr), pos(ixd1,1,1:3)
	         call gfindngb (ngb,ngbgr,ncrg,ngcr,ngus)
	        end if
          end if

	     end do
	    end do

	   end do
	
	   write(*,*) "grains' neighbours are generated"
	   do i=1,ngus
	    if (ngb(i) .ge. 1) then
	     write (40,108) i,ngb(i),(ngbgr(i,j), j=1,ngb(i))
	    end if
	   end do

      

!c-----------------------

	end program
!c--------------------
!c--------------------
!c--------------------
	   subroutine  gfindngb(ngb,ngbgr,ncrg,ngcr,ngus)
!c------This is to find a grain based on positions


	   implicit real*8 (a-h,o-z)

 	   dimension ngb(ngus),ngbgr(ngus,1000)

	   if (ngb(ncrg) .eq. 0) then
	    ngb(ncrg)=ngb(ncrg)+1
	    ngbgr(ncrg,1)=ngcr
	   else
	    idf=1
	    do i=1,ngb(ncrg)
	     if (ngbgr(ncrg,i) .eq. ngcr) idf=0 !already counted!
	    end do	
	    if (idf .eq. 1) then
	     ngb(ncrg)=ngb(ncrg)+1
	     ngbgr(ncrg,ngb(ncrg))=ngcr
	    end if

	   end if
	

	   return
	   end
!c----------------------
!c----------------------
!c----------------------

       subroutine  grainmaker (grain,nogly,pos,idl,ngus,nly,&
	   nsc,i1d,i2d,i3d) 
!C     this program generates grains based on the ceneter of the mass 
!c     of the grains and grain diameters. Data to be used for this analysis
!c     are imported from post processing on Synchrotron Data 



       implicit real*8 (a-h,o-z)

 	   dimension grain(ngus,nly,50),nogly(nly),pos(nsc,nly,4),&
	   cmsb(ngus,3),seed(ngus,3),neulprint(ngus),npgel(ngus)


	   CHARACTER*20 FILENAME
	   CHARACTER NUMBER(10)
	   DATA NUMBER/'0','1','2','3','4','5','6','7','8','9'/


107	   format(100(2x,f14.5))

	   write (*,*) 'h1'

	   do i=1,ngus
	    neulprint(i)=0
	   end do

	   idp=0
	   do i=1,i1d
	    write (*,*) i
	     do j=1,i2d
	   
	      do ik=1,i3d
	  
	       idp=idp+1
	    		
	       dist1=(pos(idp,idl,1)-grain(1,idl,4))**2+&
		         (pos(idp,idl,2)-grain(1,idl,5))**2+&
	             (pos(idp,idl,3)-grain(1,idl,6))**2-&
		         (grain(1,idl,43)/2.)**2
	  
	        do k=2,nogly(idl)
	         dist2=(pos(idp,idl,1)-grain(k,idl,4))**2+&
		          (pos(idp,idl,2)-grain(k,idl,5))**2+&
	              (pos(idp,idl,3)-grain(k,idl,6))**2-&
	              (grain(k,idl,43)/2.)**2
	    
		     if (dist2 .lt. dist1) then
	           pos(idp,idl,4)=grain(k,idl,1)
	           dist1=dist2

	         end if

	        end do

	      end do

	    end do
	   end do


	   write (*,*) 'h2',idp

	   do i=1,nogly(idl)
	    do j=1,3
	     cmsb(i,j)=0.d0
	    end do
	   end do

       cmsb=0.d0


	   do i=1,nogly(idl)
	    icgr=0
	    ngrid1=grain(i,idl,1)
	    do j=1,idp
	     if(pos(j,idl,4) .eq. ngrid1) then
	      icgr=icgr+1
	      cmsb(ngrid1,1)=cmsb(ngrid1,1)+pos(j,idl,1)
	      cmsb(ngrid1,2)=cmsb(ngrid1,2)+pos(j,idl,2)
	      cmsb(ngrid1,3)=cmsb(ngrid1,3)+pos(j,idl,3)
	     end if
	    end do
	    if (icgr .ne. 0 ) then
	     cmsb(ngrid1,1)=cmsb(ngrid1,1)/icgr
	     cmsb(ngrid1,2)=cmsb(ngrid1,2)/icgr
	     cmsb(ngrid1,3)=cmsb(ngrid1,3)/icgr
	    end if

	   end do

!CC--- This is to make sure removed grains have Center of Mass

       do i=1, nogly(idl)
        ntst=0
        ngrid1=grain(i,idl,1)
        do j=1, 3
          if (cmsb(ngrid1,j) .ne. 0) ntst=1
        end do
        if (ntst .eq. 0) then
         do j=1, 3
          cmsb(ngrid1,j)=grain(i,idl,3+j)
         end do
        end if
       
       end do



!c-----seed points
       seed=0.d0
	   do i=1,nogly(idl)
	    ngrid1=grain(i,idl,1)
	    seed(ngrid1,1)=2.d0*grain(i,idl,4)-cmsb(ngrid1,1)
	    seed(ngrid1,2)=2.d0*grain(i,idl,5)-cmsb(ngrid1,2)
	    seed(ngrid1,3)=2.d0*grain(i,idl,6)-cmsb(ngrid1,3)
	   end do

!c------working on posit

	   write (*,*) 'h3'
	   idpit=0
	   do i=1,i1d

	    write(*,*) i
	    do j=1,i2d

	     do ik=1,i3d
	    
		  idpit=idpit+1
          ngrid1=grain(1,idl,1)
	      pos(idpit,idl,4)=ngrid1
      
	      dist1=(pos(idpit,idl,1)-seed(ngrid1,1))**2+&
		        (pos(idpit,idl,2)-seed(ngrid1,2))**2+&
	            (pos(idpit,idl,3)-seed(ngrid1,3))**2-&
		        (grain(1,idl,43)/2.)**2
	  
	      do k=2,nogly(idl)
	       ngrid1=grain(k,idl,1)
	       dist2=(pos(idpit,idl,1)-seed(ngrid1,1))**2+&
		         (pos(idpit,idl,2)-seed(ngrid1,2))**2+&
	             (pos(idpit,idl,3)-seed(ngrid1,3))**2-&
	             (grain(k,idl,43)/2.)**2
	    
		   if (dist2 .lt. dist1) then
	        pos(idpit,idl,4)=grain(k,idl,1)
	        dist1=dist2

	       end if
	      end do

	     end do
	    end do
	   end do
	   write (*,*) 'h4', idpit

!c	do i=1,nogr
!c	  
!c	  nop=0
!c	  if (grain (i,9) .eq. 2) then
!c	
!c	    do j=1,nogr
!c	      if (grain(j,10) .eq. i) then
!cc	        nop=nop+1
!c			twin (i,nop+1)=j
!c	      end if
!c	    end do
!c	    twin(i,1)=nop
!c
!c	  end if
!c
!c	  write(3,107) float(i),grain(i,9),(twin(i,j),j=1,nop)
!c	end do


!c-----------------

!c	do i=1,idp
!c	  iang=pos(i,3)
!c	  write(2,107) (pos(i,j),j=1,2),(grain(iang,j),j=5,7)
!c	end do

!c	do i=1,nogr
!c	  do j=1,idp
!c	    if (pos(j,3) .eq. i) then
!c	      write(2,107) float(i), (pos(j,k),k=1,2),(grain(i,k),k=5,7)
!c	    end if
!c	  end do
!c
!c	end do


	   i1=mod(idl,10)
	   i2=(idl-i1)/10
	   i1=i1+1
	   i2=i2+1

	   FILENAME='GRAI0'//NUMBER(i2)//NUMBER(i1)//'.dat'	

	   OPEN(UNIT=50,FILE=FILENAME)
	
	   OPEN(UNIT=60,FILE='grain_euler.dat')
	   OPEN(UNIT=61,FILE='meshed_grains.dat')
       OPEN(UNIT=70,FILE='grain_euler2.dat')
	
       npgel=0.d0
       write (*,*) 'h5'
	   do i=1,nogly(idl)
	    ngrid1=grain(i,idl,1)
	    do j=1,idpit
	     if (pos(j,idl,4) .eq. ngrid1) then
	     write(50,107) float(ngrid1), (pos(j,idl,k),k=1,3),&
		               (grain(i,idl,k),k=7,9)
		               
		 npgel(i)=npgel(i)+1  !number of element per grains

		 if (neulprint(i) .eq. 0 ) then
	      neulprint(i)=1
	      write(60,107) float(ngrid1),grain(i,idl,1),&
	        (grain(i,idl,k),k=7,9)
          
          
	      write(70,107) float(ngrid1),grain(i,idl,1),&
	        (grain(i,idl,k),k=7,9), grain(i,idl,3),grain(i,idl,43)
          
	      end if
	                  
	     end if
	    end do

	   end do

	   do i=1, nogly(idl)
	    if (npgel(i) .ne. 0) then
	     write (61, 107) grain(i,idl,1), float(npgel(i))
	    end if 
	   end do


	
	   return
	   end

!cc-------------------------------------
!c--------------------------------------
!c--------------------------------------

	   subroutine rodrig2euler(qtst,angl)

!c------angle is in degree;qtst is rodrigues vect 
!c	qtst(1:3)::rodrigues vector
!c
!c     angl(1)=Phi1; 2::phi, 3::phi2 (in degree)
 
       implicit real*8 (a-h,o-z)
       dimension qtst(3),angl(3),rd(3),eul(3),rot(3,3),&
	   phi1n(2),phi2n(2),etst(3),rot1(3,3),rot2(3,3), &
       rodrigrot(3,3)

	   pi=4.d0*atan(1.d0)
	
	   rabs=qtst(1)**2+qtst(2)**2+qtst(3)**2
	   denom=1.d0/(1.d0+rabs)
	
       

     
	
	   rot(1,1)=denom*(1-rabs+2.d0*qtst(1)**2)
	   rot(1,2)=denom*(2.d0*qtst(1)*qtst(2)-2.d0*qtst(3))
	   rot(1,3)=denom*(2.d0*qtst(1)*qtst(3)+2.d0*qtst(2))
	 
	   rot(2,1)=denom*(2.d0*qtst(2)*qtst(1)+2.d0*qtst(3))
	   rot(2,2)=denom*(1-rabs+2.d0*qtst(2)**2)
	   rot(2,3)=denom*(2.d0*qtst(2)*qtst(3)-2.d0*qtst(1))

	   rot(3,1)=denom*(2.d0*qtst(1)*qtst(3)-2.d0*qtst(2))
	   rot(3,2)=denom*(2.d0*qtst(3)*qtst(2)+2.d0*qtst(1))
	   rot(3,3)=denom*(1-rabs+2.d0*qtst(3)**2)


       !!! This is to fix rigid rotation if there is any!!!
       tetx=0.d0 
       tety=0.d0 
       tetz=0.d0
       call rotaterodrig(tetx,tety,tetz,rodrigrot)
       call mulmat (rodrigrot,rot,rot1)
       
       
	   do i=1,3
	    do j=1,3
	  !   rot1(i,j)=rot(i,j)
	    end do
	   end do


	   do i=1,3
	    do j=1,3
	     rot(i,j)=rot1(j,i)
	    end do
	   end do


	   eul(2)=acos(rot(3,3))

	   if (eul(2) .le. 0.01) then
	    phi1n(1)=atan(rot(1,2)/rot(1,1))/2.d0
	    phi2n(1)=-phi1n(1)
	
	   else

	    phi2n(1)=atan((rot(1,3)/sin(eul(2)))/&
	                (rot(2,3)/sin(eul(2))))

	    phi1n(1)=atan((rot(3,1)/sin(eul(2)))/&
	                (-rot(3,2)/sin(eul(2))))
	   end if

	
	   if (phi2n(1) .ge.  0.) then
	    phi2n(2)=phi2n(1)+pi
	   else
	    phi2n(2)=phi2n(1)+pi
	    phi2n(1)=phi2n(1)+2.d0*pi
	   end if

	
	   if (phi1n(1) .ge.  0.) then
	    phi1n(2)=phi1n(1)+pi
	   else
	    phi1n(2)=phi1n(1)+pi
	    phi1n(1)=phi1n(1)+2.d0*pi
	   end if
		  
	   etst(2)=eul(2)
	   ii=1
	   jj=1
	   itest=0

	   do i=1,2

	    do j=1,2
	     etst(1)=phi1n(i)
	     etst(3)=phi2n(j)

	     call euler(etst*180.d0/pi,rot1)

!c---------It's NOT required as bothe are given in the same sense. 
!c				do ii1=1,3
!c					do jj1=1,3
!c						rot2(ii1,jj1)=rot1(ii1,jj1)
!c					end do
!c				end do
!c
!c				do ii1=1,3
!c					do jj1=1,3
!c						rot1(ii1,jj1)=rot2(jj1,ii1)
!c					end do
!c				end do

!c---------------------------
	     b=0.d0

	     do k=1,3
	      do l=1,3
	       if(abs(rot1(k,l)-rot(k,l)) .ge. 1e-7) b=1.d0
	      end do
	     end do

		 if (b .eq. 0.) then
	      itest=1
	      ii=i
	      jj=j
		 end if
	    end do
	   end do

	   if (itest .eq. 0) & 
	    write(*,*) "Euler angles are not Exact!!!"

	   eul(1)=phi1n(ii)*180.d0/pi
	   eul(3)=phi2n(jj)*180.d0/pi
	   eul(2)=eul(2)*180.d0/pi
	 
	   angl(1)=eul(1)
	   angl(2)=eul(2)
	   angl(3)=eul(3)






	   return
	   end


!c--------------------
!c--------------------
!c--------------------

       subroutine euler(prop,rotate)

!c=-----this is to go from Global to LOCAL (NO TRNASPOSE)
       implicit real*8 (a-h,o-z)
       dimension prop(3), rotate(3,3), rot(3,3)
	
!c      transforming the global coordinate to the local!
!c
!c
	   pi=4.d0*atan(1.d0)
	   phi1=prop(1)*pi/180.d0
	   phi=prop(2)*pi/180.d0
	   phi2=prop(3)*pi/180.d0
	   rotate(1,1)=cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(phi)
	   rotate(1,2)=sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(phi)
	   rotate(1,3)=sin(phi2)*sin(phi)
	   rotate(2,1)=-cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(phi)
	   rotate(2,2)=-sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(phi)
	   rotate(2,3)=cos(phi2)*sin(phi)
	   rotate(3,1)=sin(phi1)*sin(phi)
	   rotate(3,2)=-cos(phi1)*sin(phi)
	   rotate(3,3)=cos(phi)



	   return
	   end

!c------------------------
!c------------------------
!c------------------------

       subroutine rotaterodrig(tetx1,tety1,tetz1,rotmat)
!!!! Angles are in degree

       implicit real*8 (a-h,o-z)
       dimension rotmat(3,3), rotmatx(3,3), rotmaty(3,3), rotmatz(3,3), &
       rotry(3,3)
       
       pi=4.d0*atan(1.d0)
       tetx=tetx1*pi/180.d0
       tety=tety1*pi/180.d0
       tetz=tetz1*pi/180.d0
       
       rotmatx=0.d0
       rotmaty=0.d0
       rotmatz=0.d0
       
       rotmatx(1,1)=1.d0
       rotmatx(2,2)=cos(tetx)
       rotmatx(2,3)=-sin(tetx)
       rotmatx(3,2)=sin(tetx)
       rotmatx(3,3)=cos(tetx)
       
       
       rotmaty(2,2)=1.d0
       rotmaty(1,1)=cos(tety)
       rotmaty(1,3)=sin(tety)
       rotmaty(3,1)=-sin(tety)
       rotmaty(3,3)=cos(tety)
       
       
       rotmatz(3,3)=1.d0
       rotmatz(1,1)=cos(tetz)
       rotmatz(1,2)=-sin(tetz)
       rotmatz(2,1)=sin(tetz)
       rotmatz(2,2)=cos(tetz)
       
       call mulmat (rotmatx,rotmaty,rotry)
       call mulmat (rotry,rotmatz,rotmat)
       


	   return
	   end

!c------------------------
!c------------------------
!c------------------------
!c-----------------------------------------------------------
!c-----------------------------------------------------------
!c-------------------------------multiplication of two matrix
!c-----------------------------------------------------------
      subroutine mulmat (a,b,c)
      
      implicit real*8 (a-h,o-z)
      dimension a(3,3),b(3,3),c(3,3)
      
      do i=1,3
          do j=1,3
              test=0.d0
              do k=1,3
                  test=test+a(i,k)*b(k,j)
              end do
              c(i,j)=test
          end do
      end do
      
      return
      end
      
      !c----------------------------------------------------------------------
!c-----------------------------------------------------------
!c-----------------------------------------------------------
!c-----------------------------------------------------------



