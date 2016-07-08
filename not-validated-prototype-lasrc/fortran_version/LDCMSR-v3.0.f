        program ldcmsursubset
	
	implicit none
        integer(2), allocatable :: band(:,:)
        integer(2), allocatable :: aerob1(:,:)
        integer(2), allocatable :: aerob2(:,:)
        integer(2), allocatable :: aerob4(:,:)
        integer(2), allocatable :: aerob5(:,:)
        integer(2), allocatable :: aerob7(:,:)
        integer(2), allocatable :: pband(:,:)
        integer(2), allocatable :: sband(:,:,:)
	integer(2), allocatable :: tband(:,:,:)
	integer(2), allocatable :: opband(:,:)
	integer(2), allocatable :: prband(:,:)
	integer(2), allocatable :: pgband(:,:)
	integer(2), allocatable :: pbband(:,:)
	integer(2), allocatable :: oband(:,:)
	integer(2), allocatable :: tpband(:,:)
	integer(2), allocatable :: aotband(:,:)
	real tmpb10,tmpb11
	real, allocatable :: tlat(:,:)
	real, allocatable :: tlon(:,:)
	real, allocatable :: toz(:,:)
	real, allocatable :: twv(:,:)
	real, allocatable :: twvi(:,:)
	real, allocatable :: tozi(:,:)
	real, allocatable :: txcmg(:,:)
	real, allocatable :: tycmg(:,:)
	real, allocatable :: tp(:,:)
	real, allocatable :: taero(:,:)
	real, allocatable :: taero3(:,:)
	real, allocatable :: tresi(:,:)
	real, allocatable :: slprb1(:,:)
	real, allocatable :: intrb1(:,:)
	real, allocatable :: slprb2(:,:)
	real, allocatable :: intrb2(:,:)
	real, allocatable :: slprb7(:,:)
	real, allocatable :: intrb7(:,:)
	integer (2), allocatable :: tratiob1(:,:)
	integer (2), allocatable :: tratiob2(:,:)
	integer (2), allocatable :: tratiob7(:,:)
	integer (2), allocatable :: tnit(:,:)
	integer(2), allocatable :: wv(:,:)
	BYTE, allocatable :: oz(:,:)
	BYTE, allocatable :: cloud(:,:)
	BYTE, allocatable :: ipflag(:,:)
	integer(2), allocatable :: dem(:,:)
	integer(2), allocatable :: ratiob1(:,:)
	integer(2), allocatable :: ratiob2(:,:)
	integer(2), allocatable :: ratiob7(:,:)
	integer(2), allocatable :: intratiob1(:,:)
	integer(2), allocatable :: intratiob2(:,:)
	integer(2), allocatable :: intratiob7(:,:)
	integer(2), allocatable :: slpratiob1(:,:)
	integer(2), allocatable :: slpratiob2(:,:)
	integer(2), allocatable :: slpratiob7(:,:)
	integer(2), allocatable :: andwi(:,:)
	integer(2), allocatable :: sndwi(:,:)
	real th1,th2
	real*8 mall
	integer nit
	integer jj
	
	character(100) filename(12),fname,filenameanc,filenamehdf
	
	character(2) suffix(12)
	integer*8 padding
	integer sdsind(12)
	data sdsind /0,1,2,3,4,5,6,-1,7,8,9,10/
	data suffix /"01","02","03","04","05","06","07","08","09","10","11","QA"/
	integer ii,nr,nc,ib,als,i,j,nrp,ncp,ierr
	integer nrcmg,nccmg,icmg,jcmg
	real u,v
	real xcmg,ycmg
	real xndwi
	real xts,xfs,dsol,cpi
	real scalefactor,offset
	integer ihdf
	integer iband
! INPUT PARAMETER
         character*6 adate
	 integer imonth,iday,iyear        	
	 integer sfsattr
! BEGIN OF HDF PARAMETER BLOCK        
	 character*4 sdsname
         integer sfstart, sfselect, sfrdata, sfendacc, sfend,set_qamap,set_proj
	 integer sfscatt,sfwdata,sfcreate,sfginfo
         integer sd_id, sd_id2,sds_id, sds_index, status
         integer start(5), edges(5), stride(5)
         integer nstart(2), nedges(2), nstride(2)
         integer DFACC_READ,DFACC_RDWR,DFNT_CHAR8,DFACC_CREATE
	 integer DFNT_INT16,DFNT_FLOAT32,DFNT_UINT8
         parameter (DFACC_READ = 1)
         parameter (DFACC_RDWR = 3)
         parameter (DFACC_CREATE = 4)
         parameter (DFNT_INT16 = 22)
         parameter (DFNT_UINT8 = 21)
          parameter (DFNT_CHAR8 = 4)
        parameter (DFNT_FLOAT32 = 5)
	 character*80 sds_name
	 integer rank,data_type
	 integer n_attrs
	 integer dim_sizes(5)
	 integer dims(2)
         integer dim_length(5), comp_type, comp_prm(4)
! END OF HDF PARAMETER BLOCK	 

!begin block atmospheric correction variables
c   Aerosol model (string) and length of the string
       character*80 CAMOD
       integer iendarg
C Look up table for atmospheric and geometric quantities
       real tauray(16)
       real oztransa(16)
       real wvtransa(16),wvtransb(16),wvtransc(16)
       real ogtransa0(16),ogtransa1(16)
       real ogtransb0(16),ogtransb1(16)
       real ogtransc0(16),ogtransc1(16)      
       real rolutt(16,7,22,8000),nbfi(22,20),tsmax(22,20),tsmin(22,20)
       real ttv(22,20),nbfic(22,20)
       real tts(22)
       integer indts(22)
       real transt(16,7,22,22)
       real sphalbt(16,7,22),normext(16,7,22)
       real xtsmin,xtsstep,xtvmin,xtvstep,pi
       real aot550nm(22)
       data aot550nm /0.01,0.05,0.10,0.15,0.20,0.30,0.40,0.60,
     &                0.80,1.00,1.20,1.40,1.60,1.80,2.00,
     &                2.30,2.60,3.00,3.50,4.00,4.50,5.00/
       real tpres(7)
       data tpres/1050.,1013.0,900.0,800.0,700.0,600.0,500.0/
c usefull variables       
       character*6 sbandname(16)
       character*80 err_msg
C The following arguments are all names of the LUTs to look up.
       character*256 tauraynm
       character*256 gscoefnm
       character*256 anglehdf
       character*256 intrefnm
       character*256 transmnm
       character*256 spheranm
        real rotoa,xtv,xfi,raot550nm,uoz,uwv,pres
!	integer ib
        integer retval
	real roslamb
       data (sbandname(i),i=1,8)/"ldcmb1","ldcmb2","ldcmb3","ldcmb4","ldcmb5","ldcmb6",
     s  "ldcmb7","ldcmb8"/
        real tgo,roatm,ttatmg,satm,xrorayp,next
c       integer ldcmind(9)
c       data ldcmind /9,10,4,1,2,6,7,4,6/ 
!end 	block atmospheric correction variables
!block sharpening
       real r,g,b,xi,xh,xs
!end block sharpening

       real x0,y0,x,y,gsize,lat,lon
       integer row0,col0,utmzone,row,col
       real pres11,pres12,pres21,pres22
       integer uoz11,uoz21,uoz12,uoz22
       
       real erelc(8),troatm(16)
       real btgo(8),broatm(8),bttatmg(8),bsatm(8)
       integer iband1,iband3
       real raot,residual,rsurf
       real fac,xmus,corf
       integer k,l
       real xcals,xcalo,k1b10,k2b10,k1b11,k2b11,tmpf
       data xcals/3.3420E-04/
       data xcalo/0.10000/
       data k1b10/774.89/
       data k1b11/480.89/
       data k2b10/1321.08/
       data k2b11/1201.14/
       real lat1,lat2,lon1,lon2
       integer colp,rowp
       real dy,dx,ang
       integer nbval,nbclear
       real*8 anom,stemp,mclear
       real fack,facl,cldh
       integer cldhmin,cldhmax,icldh
       integer mband5,mband5k,mband5l
       real tcloud
       integer*2 cir,cld,clda,snow,clds,cldt,wat
       real cfac
c bit location of weight for cloudmask QA       
       data cld/1/
       data cir/0/
       data clda/2/
       data clds/3/
       data cldt/4/
       data wat/7/
       data cfac/6.0/
       real aaot,sresi
       real fndvi
       integer nbaot,step,hole
       real ros4,ros5
       character*512 qa  
       integer rstep,rstepc 
       integer iverbose   
       integer frln,lrln,frcl,lrcl
       integer fpln,lpln,fpcl,lpcl
       integer nci,nri
	integer xd,yd,zone,sphere
	real uplx,uply,lorx,lory,wbc,ebc,nbc,sbc
	real ros2b1
	integer supind,infind,intst,inten,x1,x2
	real rb1,rb2
	real slpr11,slpr12,slpr21,slpr22
	real intr11,intr12,intr21,intr22
	integer isuccess
	character*200 pfileout,fileout
	character*22 sstr
	integer iout

!initialisation for look up table
       iverbose=0
       xtv=0.
       xfi=0.
       xtsmin=0
       xtsstep=4.0
       xtvmin=2.84090
       xtvstep=(6.52107-2.84090)
	CAMOD="URBANCLEAN-V2.0"
	iendarg=10
	tauraynm="LDCMLUT/tauray-ldcm.ASC"
	gscoefnm="LDCMLUT/gascoef-ldcm.ASC"
	anglehdf="LDCMLUT/ANGLE_NEW.hdf"
	intrefnm="LDCMLUT/RES_LUT_V3.0-URBANCLEAN-V2.0.hdf"
	transmnm="LDCMLUT/TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII"
	spheranm="LDCMLUT/AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII"
	call readluts(CAMOD,iendarg,tauray,oztransa,wvtransa,
     s                      wvtransb,wvtransc,ogtransa0,ogtransa1,
     s                      ogtransb0,ogtransb1,ogtransc0,ogtransc1,
     s		            tsmax,tsmin,ttv,tts,nbfi,nbfic,indts,
     s                      rolutt,transt,sphalbt,normext,sbandname,err_msg,
     s                      retval,tauraynm,gscoefnm,anglehdf,intrefnm,
     s                      transmnm, spheranm)
        write(6,*) "the luts for urban clean case v2.0 have been read"
	write(6,*) "we can now perform atmospheric correction"
	
c read the DEM
	nc=7200
	nr=3600
       allocate (dem(nc,nr),stat=ierr)
       sd_id= sfstart("CMGDEM.hdf",DFACC_READ)
       start(1)=0
       start(2) = 0
       edges(1) = nc
       edges(2) = nr
       stride(1) = 1
       stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,dem)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       write(6,*) "DEM READ ", dem(2001,1001)
       
c read the RATIOFILE
	nc=7200
	nr=3600
       allocate (ratiob1(nc,nr),stat=ierr)
       allocate (ratiob2(nc,nr),stat=ierr)
       allocate (ratiob7(nc,nr),stat=ierr)
       allocate (intratiob1(nc,nr),stat=ierr)
       allocate (intratiob2(nc,nr),stat=ierr)
       allocate (intratiob7(nc,nr),stat=ierr)
       allocate (slpratiob1(nc,nr),stat=ierr)
       allocate (slpratiob2(nc,nr),stat=ierr)
       allocate (slpratiob7(nc,nr),stat=ierr)
       allocate (andwi(nc,nr),stat=ierr)
       allocate (sndwi(nc,nr),stat=ierr)
      sd_id= sfstart("ratiomapndwiexp.hdf",DFACC_READ)
      start(1)=0
       start(2) = 0
       edges(1) = nc
       edges(2) = nr
       stride(1) = 1
       stride(2) = 1

       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 1
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 4
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
      
       sds_index = 14
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,sndwi)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 18
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 19
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob1)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status








       sds_index = 15
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 16
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob2)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       
       
       
      
       sds_index = 27
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,slpratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 28
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,intratiob7)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
       sds_index = 6
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,andwi)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
c cdo not ompute ratio based on averaged ndwi never!  
c       do i=1,nc
c       do j=1,nr
c       ratiob1(i,j)=int(dble(andwi(i,j)*dble(slpratiob1(i,j)/1000.))+intratiob1(i,j))  
c       ratiob2(i,j)=int(dble(andwi(i,j)*dble(slpratiob2(i,j)/1000.))+intratiob2(i,j))  
c       ratiob7(i,j)=int(dble(andwi(i,j)*dble(slpratiob7(i,j)/1000.))+intratiob7(i,j))
c       enddo
c       enddo  
       
       write(6,*) "RATIO READ ", ratiob1(2001,1001),ratiob2(2001,1001),ratiob7(2001,1001)
       write(6,*) "RATIO READ ", slpratiob1(2001,1001),slpratiob2(2001,1001),slpratiob7(2001,1001)
       write(6,*) "RATIO READ ", intratiob1(2001,1001),intratiob2(2001,1001),intratiob7(2001,1001)
       write(6,*) "RATIO READ ", andwi(2001,1001)
c       stop
	


C Read ozone and water vapor

	pi=atan(1.)*4.
	cpi=pi
	read(5,*) ihdf
	if (ihdf.ne.0) then
	read(5,'(A100)') filenamehdf
	read(5,'(A100)') filenameanc
	read(5,*) xts,xfs
	read(5,*) utmzone,row0,col0,y0,x0
	frln=1
	frcl=1
	write(fileout,'(A200)') "correcteddata.hdf"
c
	else
c the input format is tiff
	do i=1,12
	read(5,'(A100)') filename(i)
	enddo
	read(5,'(A100)') filenameanc
	read(5,*) nr,nc,nrp,ncp
	read(5,*) xts,xfs
	read(5,*) utmzone,row0,col0,y0,x0
	read(5,'(A200)') pfileout
	read(5,*,end=19) frln,lrln,frcl,lrcl
 19     continue
        if (frln.eq.0) then
	write(fileout,'(A200)') pfileout
	frln=1
	lrln=nr
	frcl=1
	lrcl=nc
	iout=0
	else
	ii=index(pfileout," ")-5
	write(sstr,'(A2,4(A1,I4.4))') "ss","-",frln,"-",lrln,"-",frcl,"-",lrcl
	write(fileout,*) pfileout(1:ii),sstr,".hdf"
	iout=2
	endif
c	read(5,*) fpln,lpln,fpcl,lpcl
	endif
c	ii=index(pfileout," ")-1
c	jj=index(pfileout," ",.true.)
c	write(6,*) fileout(jj:ii),ii,jj,fileout
c	stop
	
	gsize=30.
	raot550nm=0.06
        fac=pi/180.
	xmus=cos(xts*fac) 




	nrcmg=3600
	nccmg=7200
       allocate (oz(nccmg,nrcmg),stat=ierr)
       allocate (wv(nccmg,nrcmg),stat=ierr)
       ii=index(filenameanc," ")-1
        fname=filenameanc(1:ii)
       sd_id= sfstart(fname,DFACC_READ)
      start(1)=0
       start(2) = 0
       edges(1) = nccmg
       edges(2) = nrcmg
       stride(1) = 1
       stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,oz)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,wv)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
c close HDF file
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
       
       write(6,*) "Ozone Water vapor read"

c	
c	allocate(band(nr,nc),STAT=als)
c	allocate(pband(nrp,ncp),STAT=als)
c	allocate(prband(nrp,ncp),STAT=als)
c	allocate(pgband(nrp,ncp),STAT=als)
c	allocate(pbband(nrp,ncp),STAT=als)
c	allocate(opband(ncp,nrp),STAT=als)
c	allocate(sband(12,nr,nc),STAT=als)
c	allocate(tband(12,nr,nc),STAT=als)
c	allocate(oband(nc,nr),STAT=als)
c	allocate(aotband(nc,nr),STAT=als)
! Getting parameter for atmospheric correction	
!        write(6,*) " aot550nm,pressure [Millibars] ,uoz [cm.atm],uwv [g/cm2]"
C update to get the parameter of the scene center
        raot550nm=0.12
	pres=1013.0
	uoz=0.30
	uwv=0.5
	xtv=0.
	xfi=0.
 	write(6,*) raot550nm,pres,uoz,uwv
	write(6,*) xts,xtv,xfi
       
c open HDF subset
        if  (ihdf.ne.0) then
        ii=index(filenamehdf," ")-1
        fname=filenamehdf(1:ii)
        sd_id= sfstart(fname,DFACC_READ)
        write(6,*) "sd_id",sd_id
        sds_index = 0
        sds_id    = sfselect(sd_id, sds_index)
        write(6,*) "sds_id", sds_id
        status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,n_attrs)
       write(6,*) "sdsname ",sds_name
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       write(6,*) dim_sizes(1),dim_sizes(2)
       nc= dim_sizes(1)
       nr=dim_sizes(2)
	frln=1
	lrln=nr
	frcl=1
	lrcl=nc
	row0=row0-nr/2
	col0=col0-nc/2
       endif 
       
       
       nci=nc
       nri=nr
       allocate(band(nci,nri),STAT=als)
       write(6,*) "allocate status ",als
c change nc,nr to subset size
       nr=lrln-frln+1
       nc=lrcl-frcl+1
       rstep=nc/50
       allocate(aerob1(nc,nr),STAT=als)
       write(6,*) "allocate status ",als
       allocate(aerob2(nc,nr),STAT=als)
       write(6,*) "allocate status ",als
       allocate(aerob4(nc,nr),STAT=als)
       write(6,*) "allocate status ",als
       allocate(aerob5(nc,nr),STAT=als)
       write(6,*) "allocate status ",als
       allocate(aerob7(nc,nr),STAT=als)
       write(6,*) "allocate status ",als
       allocate(sband(12,nc,nr),STAT=als)
       allocate(tband(12,nc,nr),STAT=als)
       allocate(oband(nc,nr),STAT=als)
       allocate(tpband(nc,nr),STAT=als)
       allocate(aotband(nc,nr),STAT=als)
       allocate(tlat(nc,nr),STAT=als)
       allocate(tlon(nc,nr),STAT=als)
       allocate(txcmg(nc,nr),STAT=als)
       allocate(tycmg(nc,nr),STAT=als)
       allocate(twv(nc,nr),STAT=als)
       allocate(twvi(nc,nr),STAT=als)
       allocate(tozi(nc,nr),STAT=als)
       allocate(tp(nc,nr),STAT=als)
       allocate(taero(nc,nr),STAT=als)
       allocate(taero3(nc,nr),STAT=als)
       allocate(tresi(nc,nr),STAT=als)
       allocate(tratiob1(nc,nr),STAT=als)
       allocate(tratiob2(nc,nr),STAT=als)
       allocate(tratiob7(nc,nr),STAT=als)
       allocate(tnit(nc,nr),STAT=als)
       allocate(cloud(nc,nr),STAT=als)
       allocate(ipflag(nc,nr),STAT=als)
       allocate(slprb1(nc,nr),STAT=als)
       allocate(intrb1(nc,nr),STAT=als)
       allocate(slprb2(nc,nr),STAT=als)
       allocate(intrb2(nc,nr),STAT=als)
       allocate(slprb7(nc,nr),STAT=als)
       allocate(intrb7(nc,nr),STAT=als)
c initialize cloud mask
      do i=1,nc
      do j=1,nr
      cloud(i,j)=0
      enddo
      enddo
       
       if (ihdf.ne.0) then     
       start(1)=0
       start(2) = 0
       edges(1) = nci
       edges(2) = nri
       stride(1) = 1
       stride(2) = 1

c read band 12 first
        ib=12
	sds_index=sdsind(ib)
        sds_id =sfselect(sd_id, sds_index)
        write(6,*) "sds_id", sds_id
        status = sfrdata(sds_id, start, stride, edges,band)
         write(6,*) "status", status
         status=sfendacc(sds_id)
         write(6,*) "status sfendacc ",status
	 else
	 ib=12
	 ii=index(filename(ib)," ")-1
	 write(6,*) "reading ",filename(ib)(1:ii)
	  open(1,file=filename(ib)(1:ii),form='UNFORMATTED',action='READ',access='DIRECT',recl=2*nci*nri+8)
	  read(1,rec=1) padding,((band(i,j),i=1,nci),j=1,nri)
	  close(1)
	 endif
	 
         do i=frcl,lrcl
	 do j=frln,lrln
	 k=i-frcl+1
	 l=j-frln+1
	 sband(ib,k,l)=band(i,j)
c	 if ((.NOT.btest(band(i,j),4)).and.(btest(band(i,j),5))) then
c	 write(6,*) " could be water "
c	 cloud(i,j)=-128
c	 endif
	 enddo
	 enddo
	 
	 
c using scene center to compute atmospheric parameter
	 row=frln+nr/2+row0
	 col=frcl+nc/2+col0
c 
     	 
	 call utmtodeg(utmzone,row,col,x0,y0,gsize,lat,lon)
	 ycmg=(89.975-lat)/0.05+1.
	 xcmg=(179.975+lon)/0.05+1
	 icmg=int(ycmg+0.5)
	 jcmg=int(xcmg+0.5)
         if (wv(jcmg,icmg).ne.0) then
	 uwv=wv(jcmg,icmg)/200.
	 else
	 uwv=0.5
	 endif
	 
         if (oz(jcmg,icmg).ne.0) then
	 uoz=oz(jcmg,icmg)/400.
	 else
	 uoz=0.3
	 endif
	 
         if (dem(jcmg,icmg).ne.-9999) then
	 pres=1013.*exp(-dem(jcmg,icmg)/8500.)
	 else
	 pres=1013.
	 endif
	 raot550nm=0.05
        	
	do ib=1,11
	if (ib.lt.9) then
	iband=ib
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
          
         endif
!computing parameter for atmospheric correction	
	if (ib.ne.8) then
	   if (ihdf.ne.0) then
	   sds_index=sdsind(ib)
           sds_id =sfselect(sd_id, sds_index)
           write(6,*) "sds_id", sds_id
           status = sfrdata(sds_id, start, stride, edges,tpband)
           write(6,*) "status", status
           status=sfendacc(sds_id)
           write(6,*) "status sfendacc ",status
           else
 	  ii=index(filename(ib)," ")-1
	  write(6,*) "reading ",filename(ib)(1:ii)
	  open(1,file=filename(ib)(1:ii),form='UNFORMATTED',action='READ',access='DIRECT',recl=2*nci*nri+8)
	  read(1,rec=1) padding,((band(i,j),i=1,nci),j=1,nri)
	  close(1)
c calibrate band 1 to 9 except 8
          if (ib.le.9) then
             do i=frcl,lrcl
	     do j=frln,lrln
 	     if (band(i,j).ge.0) then
             rotoa=(band(i,j)*2.0000E-05)-0.1
	     else
	     rotoa=((65536+band(i,j))*2.0000E-05)-0.1
	      endif
     	     k=i-frcl+1
	     l=j-frln+1
	     tpband(k,l)=int(rotoa*10000./xmus)
	     enddo
	     enddo
	     endif  
          endif
	  
	  
c call the atmospheric correction
             if (ib.lt.9) then	  
	     do i=1,nc
	     do j=1,nr
	     if (sband(12,i,j).ne.1) then
             rotoa=tpband(i,j)/10000.
	     if (ib.eq.1) aerob1(i,j)=tpband(i,j)
	     if (ib.eq.2) aerob2(i,j)=tpband(i,j)
	     if (ib.eq.4) aerob4(i,j)=tpband(i,j)
	     if (ib.eq.5) aerob5(i,j)=tpband(i,j)
	     if (ib.eq.7) aerob7(i,j)=tpband(i,j)
             roslamb=rotoa/tgo
 	     roslamb=roslamb-roatm
 	     roslamb=roslamb/ttatmg
 	     roslamb=roslamb/(1.+satm*roslamb)
	     sband(ib,i,j)=int(roslamb*10000.)
	     btgo(ib)=tgo
	     broatm(ib)=roatm
	     bttatmg(ib)=ttatmg
	     bsatm(ib)=satm
	     else
	     sband(ib,i,j)=-1000
	     endif
	     enddo
	     enddo
	     else
	     if (ihdf.ne.0) then
             do i=frcl,lrcl
	     do j=frln,lrln
     	     k=i-frcl+1
	     l=j-frln+1
	     sband(ib,k,l)=tpband(i,j)
	     enddo
	     enddo
	     else
             do i=frcl,lrcl
	     do j=frln,lrln
     	     k=i-frcl+1
	     l=j-frln+1
	     sband(ib,k,l)=band(i,j)
	     enddo
	     enddo
	     endif
	   endif
	   
	   if (ib.eq.9) then
	   do i=1,nc
	   do j=1,nr
	   sband(ib,i,j)=tpband(i,j)
	   enddo
	   enddo
	   endif
	     
	     
	     
	else
c iband eq 8	
	write(6,*) "not reading band8 sorry"
	endif   
	
c
	enddo 
	
c calibrate thermal bands
         do i=1,nc
	 do j=1,nr
	 if (sband(10,i,j).lt.0)  then 
	    tmpb10=65536+sband(10,i,j)
	    else
	    tmpb10=sband(10,i,j)
	    endif
	    
	 if (sband(11,i,j).lt.0)  then 
	     tmpb11=65536+sband(11,i,j)
	     else
	     tmpb11=sband(11,i,j)
	     endif
	     
	 if ((sband(12,i,j).ne.1)) then
	 tmpf=xcals*tmpb10+xcalo
	 tmpf=k2b10/log(1.+k1b10/tmpf)
	 sband(10,i,j)=int(tmpf*10.)
	 tmpf=xcals*tmpb11+xcalo
	 tmpf=k2b11/log(1.+k1b11/tmpf)
	 sband(11,i,j)=int(tmpf*10.)
	 endif
	 enddo
	 enddo
         write(6,*) "Finished calibrating Thermal"	 
c compute lat,lon,rowcmg,colcmg,ozone,water vapor,pressure
c skipping reading the aerosol product goto 111
         goto 111
         write(6,*) "reading aerosol data"
         sd_id= sfstart('aerosoldata.hdf',DFACC_READ)
	 dim_length(1)=nc
	 dim_length(2)=nr
	 comp_type=4
	 comp_prm(1)=8
	 rank=2
	 dim_sizes(1)=nc
	 dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
       sds_index = 0
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,taero3)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 1
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,taero)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 2
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,tresi)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       sds_index = 3
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges,ipflag)
       write(6,*) "status", status
       status = sfendacc(sds_id)
       write(6,*) "status sfendacc ",status
       status = sfend(sd_id)
       write(6,*) "status sfend ",status
111    continue
         do i=1,nc
c	 write(6,*) " Processing column ",i
	 if ((i-int(i/rstep)*rstep).eq.0) then
	 write(6,*) " processing collumn ",i
	 endif
	 do j=1,nr
c	 goto 11
	 if ((sband(12,i,j).eq.1)) then
	 ipflag(i,j)=5
	 goto 11
	 endif
	 if ((j.eq.673).and.(i.eq.709)) then
	 write(6,*) "J'en tiens un"
	 endif
	 
	 row=frln+(j-1)+row0
	 col=frcl+(i-1)+col0
	 call utmtodeg(utmzone,row,col,x0,y0,gsize,lat,lon)
	 tlat(i,j)=lat
	 tlon(i,j)=lon
	 ycmg=(89.975-lat)/0.05+1.
	 xcmg=(179.975+lon)/0.05+1
	 txcmg(i,j)=xcmg
	 tycmg(i,j)=ycmg
	 icmg=int(ycmg+0.5)
	 jcmg=int(xcmg+0.5)
c         twv(i,j)=wv(jcmg,icmg)
	 icmg=int(ycmg)
	 jcmg=int(xcmg)
	 u=(ycmg-icmg)
	 v=(xcmg-jcmg)
	 twvi(i,j)=wv(jcmg,icmg)*(1.-u)*(1.-v)+wv(jcmg+1,icmg)*(1.-u)*v
     s	      +wv(jcmg,icmg+1)*(1.-v)*u+wv(jcmg+1,icmg+1)*u*v
         twvi(i,j)=twvi(i,j)/100.
	 
	 
	 
	 if (oz(jcmg,icmg).lt.0) then
	 uoz11=256+oz(jcmg,icmg)
	 else
	 uoz11=oz(jcmg,icmg)
	 endif
	 
	 if (oz(jcmg+1,icmg).lt.0) then
	 uoz12=256+oz(jcmg+1,icmg)
	 else
	 uoz12=oz(jcmg+1,icmg)
	 endif
	 
	 if (oz(jcmg,icmg+1).lt.0) then
	 uoz21=256+oz(jcmg,icmg+1)
	 else
	 uoz21=oz(jcmg,icmg+1)
	 endif
	 
	 if (oz(jcmg+1,icmg+1).lt.0) then
	 uoz22=256+oz(jcmg+1,icmg+1)
	 else
	 uoz22=oz(jcmg+1,icmg+1)
	 endif
	 
	 if (uoz11.eq.0) uoz11=120
	 if (uoz12.eq.0) uoz12=120
	 if (uoz21.eq.0) uoz21=120
	 if (uoz22.eq.0) uoz22=120
	 
	  
	 
	 
	 tozi(i,j)=uoz11*(1.-u)*(1.-v)+uoz12*(1.-u)*v
     s	      +uoz21*(1.-v)*u+uoz22*u*v
     
         tozi(i,j)=tozi(i,j)/400.
	 
	 if (dem(jcmg,icmg).ne.-9999) then
	 pres11=1013.*exp(-dem(jcmg,icmg)/8500.)
	 else
	 pres11=1013.
	 cloud(i,j)=-128
c	 tresi(i,j)=-1.
	 endif
	 if (dem(jcmg+1,icmg).ne.-9999) then
	 pres12=1013.*exp(-dem(jcmg+1,icmg)/8500.)
	 else
	 pres12=1013.
	 endif
	 if (dem(jcmg,icmg+1).ne.-9999) then
	 pres21=1013.*exp(-dem(jcmg,icmg+1)/8500.)
	 else
	 pres21=1013.
	 endif
	 if (dem(jcmg+1,icmg+1).ne.-9999) then
	 pres22=1013.*exp(-dem(jcmg+1,icmg+1)/8500.)
	 else
	 pres22=1013.
	 endif
	 tp(i,j)=pres11*(1.-u)*(1.-v)+pres12*(1.-u)*v
     s	      +pres21*(1.-v)*u+pres22*u*v
c testing for water     
        if (1.eq.1) then
c inverting aerosols	
c filter cirrus 
        if (sband(9,i,j).gt.(100./(tp(i,j)/1013.))) then
	cloud(i,j)=cloud(i,j)+1
	else
c inverting aerosol	       
	do ib=1,8
	erelc(ib)=-1.
	enddo
	rb1=ratiob1(jcmg,icmg)/1000.
	rb2=ratiob2(jcmg,icmg)/1000.
	if ((abs((rb1-0.454878*rb2*rb2-0.459559*rb2)).gt.0.15).or.(ratiob7(jcmg,icmg).lt.1000)) then
	slpratiob1(jcmg,icmg)=0
	slpratiob2(jcmg,icmg)=0
	slpratiob7(jcmg,icmg)=0
	intratiob1(jcmg,icmg)=327
	intratiob2(jcmg,icmg)=482
	intratiob7(jcmg,icmg)=2000
	else
	if (sndwi(jcmg,icmg).lt.200) then
	slpratiob1(jcmg,icmg)=0
	slpratiob2(jcmg,icmg)=0
	slpratiob7(jcmg,icmg)=0
	intratiob1(jcmg,icmg)=ratiob1(jcmg,icmg)
	intratiob2(jcmg,icmg)=ratiob2(jcmg,icmg)
	intratiob7(jcmg,icmg)=ratiob7(jcmg,icmg)
	endif
	endif
	
	rb1=ratiob1(jcmg+1,icmg)/1000.
	rb2=ratiob2(jcmg+1,icmg)/1000.
	if ((abs((rb1-0.454878*rb2*rb2-0.459559*rb2)).gt.0.15).or.(ratiob7(jcmg+1,icmg).lt.1000)) then
	slpratiob1(jcmg+1,icmg)=0
	slpratiob2(jcmg+1,icmg)=0
	slpratiob7(jcmg+1,icmg)=0
	intratiob1(jcmg+1,icmg)=327
	intratiob2(jcmg+1,icmg)=482
	intratiob7(jcmg+1,icmg)=2000
	else
	if (sndwi(jcmg+1,icmg).lt.200) then
	slpratiob1(jcmg+1,icmg)=0
	slpratiob2(jcmg+1,icmg)=0
	slpratiob7(jcmg+1,icmg)=0
	intratiob1(jcmg+1,icmg)=ratiob1(jcmg+1,icmg)
	intratiob2(jcmg+1,icmg)=ratiob2(jcmg+1,icmg)
	intratiob7(jcmg+1,icmg)=ratiob7(jcmg+1,icmg)
	endif
	endif
	
	rb1=ratiob1(jcmg,icmg+1)/1000.
	rb2=ratiob2(jcmg,icmg+1)/1000.
	if ((abs((rb1-0.454878*rb2*rb2-0.459559*rb2)).gt.0.15).or.(ratiob7(jcmg,icmg+1).lt.1000)) then
	slpratiob1(jcmg,icmg+1)=0
	slpratiob2(jcmg,icmg+1)=0
	slpratiob7(jcmg,icmg+1)=0
	intratiob1(jcmg,icmg+1)=327
	intratiob2(jcmg,icmg+1)=482
	intratiob7(jcmg,icmg+1)=2000
	else
	if (sndwi(jcmg,icmg+1).lt.200) then
	slpratiob1(jcmg,icmg+1)=0
	slpratiob2(jcmg,icmg+1)=0
	slpratiob7(jcmg,icmg+1)=0
	intratiob1(jcmg,icmg+1)=ratiob1(jcmg,icmg+1)
	intratiob2(jcmg,icmg+1)=ratiob2(jcmg,icmg+1)
	intratiob7(jcmg,icmg+1)=ratiob7(jcmg,icmg+1)
	endif
	endif

	rb1=ratiob1(jcmg+1,icmg+1)/1000.
	rb2=ratiob2(jcmg+1,icmg+1)/1000.
	if ((abs((rb1-0.454878*rb2*rb2-0.459559*rb2)).gt.0.15).or.(ratiob7(jcmg+1,icmg+1).lt.1000)) then
	slpratiob1(jcmg+1,icmg+1)=0
	slpratiob2(jcmg+1,icmg+1)=0
	slpratiob7(jcmg+1,icmg+1)=0
	intratiob1(jcmg+1,icmg+1)=327
	intratiob2(jcmg+1,icmg+1)=482
	intratiob7(jcmg+1,icmg+1)=2000
	else
	if (sndwi(jcmg+1,icmg+1).lt.200) then
	slpratiob1(jcmg+1,icmg+1)=0
	slpratiob2(jcmg+1,icmg+1)=0
	slpratiob7(jcmg+1,icmg+1)=0
	intratiob1(jcmg+1,icmg+1)=ratiob1(jcmg+1,icmg+1)
	intratiob2(jcmg+1,icmg+1)=ratiob2(jcmg+1,icmg+1)
	intratiob7(jcmg+1,icmg+1)=ratiob7(jcmg+1,icmg+1)
	endif
	endif

	 slpr11=slpratiob2(jcmg,icmg)/1000.
	 intr11=intratiob2(jcmg,icmg)/1000.
	 slpr12=slpratiob2(jcmg+1,icmg)/1000.
	 intr12=intratiob2(jcmg+1,icmg)/1000.
	 slpr21=slpratiob2(jcmg,icmg+1)/1000.
	 intr21=intratiob2(jcmg,icmg+1)/1000.
	 slpr22=slpratiob2(jcmg+1,icmg+1)/1000.
	 intr22=intratiob2(jcmg+1,icmg+1)/1000.
	 slprb2(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s	      +slpr21*(1.-v)*u+slpr22*u*v
	 intrb2(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s	      +intr21*(1.-v)*u+intr22*u*v

	 slpr11=slpratiob1(jcmg,icmg)/1000.
	 intr11=intratiob1(jcmg,icmg)/1000.
	 slpr12=slpratiob1(jcmg+1,icmg)/1000.
	 intr12=intratiob1(jcmg+1,icmg)/1000.
	 slpr21=slpratiob1(jcmg,icmg+1)/1000.
	 intr21=intratiob1(jcmg,icmg+1)/1000.
	 slpr22=slpratiob1(jcmg+1,icmg+1)/1000.
	 intr22=intratiob1(jcmg+1,icmg+1)/1000.
	 slprb1(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s	      +slpr21*(1.-v)*u+slpr22*u*v
	 intrb1(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s	      +intr21*(1.-v)*u+intr22*u*v

	 slpr11=slpratiob7(jcmg,icmg)/1000.
	 intr11=intratiob7(jcmg,icmg)/1000.
	 slpr12=slpratiob7(jcmg+1,icmg)/1000.
	 intr12=intratiob7(jcmg+1,icmg)/1000.
	 slpr21=slpratiob7(jcmg,icmg+1)/1000.
	 intr21=intratiob7(jcmg,icmg+1)/1000.
	 slpr22=slpratiob7(jcmg+1,icmg+1)/1000.
	 intr22=intratiob7(jcmg+1,icmg+1)/1000.
	 slprb7(i,j)=slpr11*(1.-u)*(1.-v)+slpr12*(1.-u)*v
     s	      +slpr21*(1.-v)*u+slpr22*u*v
	 intrb7(i,j)=intr11*(1.-u)*(1.-v)+intr12*(1.-u)*v
     s	      +intr21*(1.-v)*u+intr22*u*v



        xndwi=dble(sband(5,i,j)-dble(sband(7,i,j)/2.))/dble(sband(5,i,j)+dble(sband(7,i,j)/2.))
	th1=((andwi(jcmg,icmg)+2.*sndwi(jcmg,icmg))/1000.)
	th2=((andwi(jcmg,icmg)-2.*sndwi(jcmg,icmg))/1000.)
	if (xndwi.gt.th1) xndwi=th1
	if (xndwi.lt.th2) xndwi=th2
        erelc(1)=(xndwi*dble(slprb1(i,j))+intrb1(i,j)) 
        erelc(2)=(xndwi*dble(slprb2(i,j))+intrb2(i,j)) 
        erelc(7)=(xndwi*dble(slprb7(i,j))+intrb7(i,j))
 	erelc(4)=1.
	tratiob1(i,j)=int(erelc(1)*1000.)
	tratiob2(i,j)=int(erelc(2)*1000.)
	tratiob7(i,j)=int(erelc(7)*1000.)
	
	
	
	
	troatm(1)=aerob1(i,j)/10000.
 	troatm(2)=aerob2(i,j)/10000.
	troatm(4)=aerob4(i,j)/10000.
	troatm(7)=aerob7(i,j)/10000.
	iband1=4
	iband3=2 
	pres=tp(i,j)
	uoz=tozi(i,j)
	uwv=twvi(i,j)
	if ((i.eq.184).and.(j.eq.151)) then
	write(6,*) "aeroband1 ",aerob1(i,j)
       write(6,*) i,j,xts,xtv,xfi,pres,uoz,uwv,erelc(1),erelc(2),erelc(7),
     s	troatm(1),troatm(2),troatm(4),troatm(7),iband1,iband3,tratiob1(i,j)
        iverbose=1
	else
	iverbose=0
        endif
c skipping the aerosol inversion goto 11
c       goto 11

c       if (btest(cloud(i,j),wat)) then
c       fndvi=dble(sband(5,i,j)-sband(4,i,j))/dble(sband(5,i,j)+sband(4,i,j))
c       if (fndvi.lt.0.1) then
c       taero(i,j)=0.05
c       tresi(i,j)=0.0
c       ipflag(i,j)=4
c       goto 11
c       endif
c       endif
       call subaeroret(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iverbose,raot,residual,nit,next)
     
c         write(6,*) "aero retrieved " ,i,j,raot,residual
c	 endif
         
         tnit(i,j)=nit
	 taero3(i,j)=next*raot
         corf=raot/xmus
c         if (residual.lt.((0.015+0.005*corf)*aerob7(i,j)/300.)) then
         if (residual.lt.((0.010+0.005*corf)*1.)) then
c test if band5 makes sense
	 
	iband=5
	rotoa=aerob5(i,j)/10000.
	raot550nm=raot
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros5=roslamb
	iband=4
	rotoa=aerob4(i,j)/10000.
	raot550nm=raot
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
         ros4=roslamb
         if ((ros5.gt.0.1).and.((ros5-ros4)/(ros4+ros5).gt.0)) then
	 taero(i,j)=raot
	 tresi(i,j)=residual
	 else
	 taero(i,j)=raot
	 tresi(i,j)=residual
	 ipflag(i,j)=2
	 endif
	 
	 else
	 taero(i,j)=raot
	 tresi(i,j)=residual
	 ipflag(i,j)=3
	 endif
c endif of else cirrus.	 
	 endif
c endif not water pixel	 
	 endif
c here to skip fill value

 11      continue	 
	 enddo
	 enddo
c saved the aerosol product,residual and ipflag for future debugging	
         goto 777 
         write(6,*) "saving aerosol data"
         sd_id= sfstart('aerosoldata.hdf',DFACC_CREATE)
	 dim_length(1)=nc
	 dim_length(2)=nr
	 comp_type=4
	 comp_prm(1)=8
	 rank=2
	 dim_sizes(1)=nc
	 dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
	 write(sds_name,'(A11)') "taero-band1"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,taero3)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	 write(sds_name,'(A5)') "taero"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,taero)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	 write(sds_name,'(A5)') "tresi"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tresi)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
         write(sds_name,'(A6)') "IPFLAG"
	 sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,ipflag)
	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         status = sfend(sd_id) 	 
	 write(6,*) "status sfend ",status
 777     continue	 
	 
	 
	 
c refine the cloud mask 	
c compute the average temperature of the clear non water non filled data
           nbval=0
	   nbclear=0  
	   mclear=0.
	   mall=0.
	  
           do i=1,nc
	   do j=1,nr
	   if (sband(12,i,j).ne.1) then
	   nbval=nbval+1
	   mall=mall+(sband(10,i,j)/10.)
	   if ((.not.btest(cloud(i,j),cir)).and.(sband(5,i,j).gt.300)) then
	   anom=sband(2,i,j)-sband(4,i,j)/2.
	   if (anom.lt.300) then
	   if ((sband(10,i,j).lt.2500)) then
	   write(6,*) "J'en tiens un autre ",i,j
	   endif
	   nbclear=nbclear+1
	   mclear=mclear+(sband(10,i,j)/10.)
	   endif
	   endif
	   endif
	   enddo
	   enddo
	   
	   if (nbclear.gt.0) then
	   mclear=mclear/nbclear
	   else
	   mclear=275.0
	   endif
	   
	   if (nbval.gt.0) then
	   mall=mall/nbval
	   endif
	   
	   write(6,*) "average clear temperature  %clear", mclear,nbclear*100./(nr*nc),mall,nbval
	   
	   

c Cloud mask 
c skipping cloud mask (goto 999)
c          goto 999
	  do i=1,nc
	  do j=1,nr
	  if (ipflag(i,j).gt.1) then
C all bad retrieval (ipflag ne 0) except for filled value (ip.eq.5) should be tested
C 
          if (((sband(2,i,j)-sband(4,i,j)/2.).gt.500).and.((sband(10,i,j)/10.).lt.(mclear-2.))) then
c snow or cloud cloud for now
          cloud(i,j)=cloud(i,j)+2
	  endif
	  endif
	  enddo
	  enddo	
	  
c	  goto 999
c set up the adjacent to something bad (snow or cloud) bit
       do i=1,nc
       do j=1,nr
       if (btest(cloud(i,j),cld).or.(btest(cloud(i,j),cir))) then
       do k=i-5,i+5
       do l=j-5,j+5
       if ((k.ge.1).and.(k.le.nc).and.(l.ge.1).and.(l.le.nr)) then
       if ((.NOT.btest(cloud(k,l),cld)).and.(.not.btest(cloud(k,l),cir)).and.(.not.btest(cloud(k,l),clda))) then
          cloud(k,l)=cloud(k,l)+4
	  endif
	endif
	enddo
	enddo
	endif
	enddo
	enddo
c cloud shadow
c compute adjustement to the North
c scene center lat,lon
	 row=frln+(nr-1)/2.
	 col=frcl+(nc-1)/2.
	 write(6,*) "row col ",row,col
	 call utmtodeg(utmzone,row,col,x0,y0,gsize,lat1,lon1)
c move 100 pixels north
         write(6,*) "LAT1 LON1",lat1,lon1
         rowp=row-100	 
	 call utmtodeg(utmzone,rowp,col,x0,y0,gsize,lat2,lon2)
         write(6,*) "LAT2 LON2",lat2,lon2
	 lon=lon1
	 lat=lat2
	 call degtoutm(utmzone,lat,lon,x0,y0,gsize,colp,rowp)
	 write(6,*) "colp, rowp",colp,rowp
	 dy=row-rowp
	 dx=colp-col
	 ang=atan(dx/dy)*180./pi
	 write(6,*) "Adjustement to true North ", ang
c  compute cloud shadow  
       write(6,*) "looking for cloud shadow"
       fack=cos(xfs*fac)*tan(xts*fac)/gsize
       facl=sin(xfs*fac)*tan(xts*fac)/gsize
       do i=1,nc
       do j=1,nr
       if (btest(cloud(i,j),cld).or.btest(cloud(i,j),cir)) then
           tcloud=sband(10,i,j)/10.
           cldh=(mclear-tcloud)*1000./cfac
	   if (cldh.lt.0.) cldh=0.
	   cldhmin=cldh-1000.
	   cldhmax=cldh+1000.
	   mband5=9999
	   if (cldhmin.lt.0) cldhmin=0.
           do icldh=cldhmin/10,cldhmax/10
	    cldh=icldh*10.
            k=j+fack*cldh
            l=i-facl*cldh
            if ((k.ge.1).and.(k.le.nr).and.(l.ge.1)
     &	    .and.(l.le.nc)) then

            if ((sband(6,l,k).lt.800).and.((sband(3,l,k)-sband(4,l,k)).lt.100)) then
            if (btest(cloud(l,k),cld).or.btest(cloud(l,k),cir)
     &       .or.btest(cloud(l,k),clds)) then
             continue
             else
c store the value of band6 as well as the l and k value
             if (sband(6,l,k).lt.mband5) then
	     mband5=sband(6,l,k)
	     mband5k=k
	     mband5l=l
	     endif    
             endif
             endif
             endif
	     enddo
	     if (mband5.lt.9999) then
	     l=mband5l
	     k=mband5k
             cloud(l,k)=cloud(l,k)+(2**clds)
	     endif
         endif
	 enddo
	 enddo
	 
c exgd the cloud shadow using residual 	 
c need to be checked carefully residual is meaningless here
       do i=1,nc
       do j=1,nr
       if (btest(cloud(i,j),clds)) then
       do k=i-6,i+6
       do l=j-6,j+6
       if ((k.ge.1).and.(k.le.nc).and.(l.ge.1).and.(l.le.nr)) then
       if (btest(cloud(k,l),cld)
     &   .or.btest(cloud(k,l),clds)) then
       continue
       else
       if (btest(cloud(k,l),cldt)) then
          continue
          else
           if (tresi(k,l).lt.0) then
           cloud(k,l)=cloud(k,l)+(2**cldt)
           endif
       endif
       endif
       endif
       enddo
       enddo
       endif
       enddo
       enddo
c update the cloud shadow       
       write(6,*) "updating cloud shadow"
       do i=1,nc
       do j=1,nr
       if (btest(cloud(i,j),cldt)) then
       cloud(i,j)=cloud(i,j)+(2**clds)
       cloud(i,j)=cloud(i,j)-(2**cldt)
       endif
       enddo
       enddo
c detect water (ipflag=6)
       do j=1,nr
       do i=1,nc
       if ((.not.btest(cloud(i,j),cir)).and.(.not.btest(cloud(i,j),cld))
     s    .and.(.not.btest(cloud(i,j),clds)) .and.(.not.btest(cloud(i,j),clda)).and.(sband(4,i,j).lt.1000)) then
        if (sband(5,i,j).lt.100) then
	fndvi=-0.01
	else
        fndvi=dble(sband(5,i,j)-sband(4,i,j))/dble(sband(5,i,j)+sband(4,i,j))
        endif
	if (fndvi.lt.0.01) then
        ipflag(i,j)=6
        endif
	endif
c detect snow/ confusion with high aerosol for now (flag as cloud)
        if (sband(4,i,j).gt.1000) then
	fndvi=(sband(4,i,j)-sband(7,i,j)/2.)/dble(sband(4,i,j)+sband(7,i,j)/2.)
	if (fndvi.gt.0.5) then
        taero(i,j)=0.05
        ipflag(i,j)=7
        cloud(i,j)=cloud(i,j)+2
        endif
	endif
	
	enddo
	enddo
c aerosol interpolation
999    continue
c begin skipping interpolation (goto 1001)
c       goto 1001

       do j=1,nr
       do i=1,nc
        if ((ipflag(i,j).gt.1).and.(.not.btest(cloud(i,j),cir)).and.(.not.btest(cloud(i,j),cld))
     s    .and.(.not.btest(cloud(i,j),wat)).and.(ipflag(i,j).lt.5)) then
c problem there is no interpolation over water      
c interpolation 
       
       k=i
       do while ((ipflag(k,j).gt.1).and.(k.gt.1)) 
       k=k-1     
       enddo
       infind=k
       k=i
       do while ((ipflag(k,j).gt.1).and.(k.lt.nc))
       k=k+1     
       enddo
       supind=k
       
       intst=infind+1
       inten=supind-1
       
       
       
1003   continue
       
       if ((ipflag(infind,j).gt.1).or.(ipflag(supind,j).gt.1))then
        goto 1002
       endif
       
      
       do k=intst,inten
       if (ipflag(k,j).lt.5) then
       taero(k,j)=taero(infind,j)+(taero(supind,j)-taero(infind,j))*(k-infind)/(supind-infind)
       ipflag(k,j)=1
       endif
       enddo
 1004  continue
       
       endif
 1002   continue
       enddo
       enddo
c

            
c redo in the reverse direction
c       goto 2005
       do i=1,nc
       do j=1,nr
        if ((ipflag(i,j).gt.1).and.(.not.btest(cloud(i,j),cir)).and.(.not.btest(cloud(i,j),cld))
     s    .and.(.not.btest(cloud(i,j),wat)).and.(ipflag(i,j).lt.5)) then
c problem there is no interpolation over water      
c interpolation 
       
       k=j
       do while ((ipflag(i,k).gt.1).and.(k.gt.1)) 
       k=k-1     
       enddo
       infind=k
       k=j
       do while ((ipflag(i,k).gt.1).and.(k.lt.nr))
       k=k+1     
       enddo
       supind=k
       
       intst=infind+1
       inten=supind-1
       
       
       
2003   continue
       
       if ((ipflag(i,infind).gt.1).or.(ipflag(i,supind).gt.1))then
        goto 2002
       endif
       
      
       do k=intst,inten
       if (ipflag(i,k).lt.5) then
       taero(i,k)=taero(i,infind)+(taero(i,supind)-taero(i,infind))*(k-infind)/(supind-infind)
       ipflag(i,k)=1
       endif
c       tresi(i,k)=0.0
       enddo
 2004  continue
       
       endif
 2002  continue
       enddo
       enddo
       
 2005  continue     
c

            




           	 
c perform the correction
 1001     continue
	  rstepc=rstep*5
          do ib=1,7
	  do i=1,nc
	 if ((i-int(i/rstepc)*rstepc).eq.0) then
	 write(6,*) " atmospheric correction collumn ",i," band ",ib
	 endif
	  do j=1,nr
c	  write(6,*) "i , j",i,j,tresi(i,j)
	  if (sband(12,i,j).ne.1) then
	  if ((.not.btest(cloud(i,j),cir)).and.(.not.btest(cloud(i,j),cld))) then
	  rsurf=sband(ib,i,j)/10000.
	  rotoa=(rsurf*bttatmg(ib)/(1.-bsatm(ib)*rsurf)+broatm(ib))*btgo(ib)
	  if (ib.eq.1) then
c	  write(6,*) "band 1", aerob1(i,j),rotoa
	  endif
	  raot550nm=taero(i,j)
	  pres=tp(i,j)
	  uwv=twvi(i,j)
	  uoz=tozi(i,j)
	  iband=ib
c	  write(6,*) xts,xtv,xfi,raot550nm,iband,pres,uoz,uwv,tauray
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
          if ((i.eq.2111).and.(j.eq.4123)) then
	  write(6,*) "roslamb ",roslamb," tauaero ",taero(i,j) 
	  endif
          if (ib.eq.1) then
	  if ((roslamb.lt.0.01).and.(raot550nm.gt.0.06)) then
c this should be refined to get better value	  
          isuccess=0
	  do while ((roslamb.lt.0))
          if ((i.eq.1001).and.(j.eq.237)) then
	  write(6,*) "roslamb ",roslamb," tauaero ",taero(i,j) 
	  endif
	  raot550nm=taero(i,j)-0.05
	  ros2b1=roslamb
	  pres=tp(i,j)
	  uwv=twvi(i,j)
	  uoz=tozi(i,j)
	  iband=ib
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
         taero(i,j)= raot550nm
	 if (roslamb.ge.0.) isuccess=1
	 enddo
c enddo of do while lt. 0	
 
	 if (isuccess.eq.0) then
	  raot550nm=taero(i,j)-0.05
	  ros2b1=roslamb
	  pres=tp(i,j)
	  uwv=twvi(i,j)
	  uoz=tozi(i,j)
	  iband=ib
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
         taero(i,j)= raot550nm
         endif	 
	 
c refined optical depth	 
         raot550nm=taero(i,j)-(roslamb-0.01)*(-0.05)/(roslamb-ros2b1)
	 taero(i,j)=raot550nm
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       rotoa,roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
	 endif 
c endif or refining optical depth
c set up aerosol QA
c          write(6,*) "SETTING UP aerosol qa BITS"
          if (abs(rsurf-roslamb).le.0.015) then
	    cloud(i,j)=cloud(i,j)+16	
	    else
	    if   (abs(rsurf-roslamb).lt.0.03) then
	          cloud(i,j)=cloud(i,j)+32
	          else
	          cloud(i,j)=cloud(i,j)+48
	          endif
	    endif
	  endif
C endif of ib=1	  
	  
	  sband(ib,i,j)=int(roslamb*10000.)
	  endif
	  endif
	  enddo
	  enddo
	  enddo

         
!Saving in an HDF file
         write(6,*) "saving corrected data"
	 ii=index(fileout," ")-1
	 write(6,*) "fileout ",fileout
         sd_id= sfstart(adjustl(fileout),DFACC_CREATE)
	 write(6,*) " filout ",sd_id
	 do ib=1,12
	 if (ib.ne.8) then
	 do i=1,nc
	 do j=1,nr
	 oband(i,j)=sband(ib,i,j)
	 enddo
	 enddo
	 if (ib.eq.2) then
         write(sds_name,'(A4,A2,A5)') "band",suffix(ib),"-blue"
	 endif
	 if (ib.eq.3) then
         write(sds_name,'(A4,A2,A6)') "band",suffix(ib),"-green"
	 endif
	 if (ib.eq.4) then
         write(sds_name,'(A4,A2,A4)') "band",suffix(ib),"-red"
	 endif
	 if ((ib.le.1).or.(ib.ge.5)) then
         write(sds_name,'(A4,A2)') "band",suffix(ib)
	 endif
	 write(6,*) "writing band ",sds_name
	 dim_length(1)=nc
	 dim_length(2)=nr
	 comp_type=4
	 comp_prm(1)=8
	 rank=2
	 dim_sizes(1)=nc
	 dim_sizes(2)=nr
         start(1) = 0
         start(2) = 0
          edges(1) = nc
         edges(2) = nr
         stride(1) = 1
         stride(2) = 1
	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,oband)
c check if QA band to write the attribute
	 write(6,*) "status sfwdata ",status
          if ((ib.ge.1).and.(ib.le.9)) then
	 write(6,*) "writing attribute band ",sds_id
	 scalefactor=10000.
	 offset=0.0
	 status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
	 status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
	endif
          if ((ib.eq.10).or.(ib.eq.11)) then
	 write(6,*) "writing attribute band ",sds_id
	 scalefactor=10.
	 offset=0.0
	 status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
	 status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
	endif
        status = sfendacc(sds_id)
	 write(6,*) "status sfendacc ",status
	 endif
	 enddo
c writing lat and lon
c            write(sds_name,'(A3)') "lat"
c   	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
c            status=sfwdata(sds_id,start,stride,edges,tlat)
c   	 write(6,*) "status sfwdata ",status
c            status = sfendacc(sds_id)
c            write(sds_name,'(A3)') "lon"
c   	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
c            status=sfwdata(sds_id,start,stride,edges,tlon)
c   	 write(6,*) "status sfwdata ",status
c            status = sfendacc(sds_id)
c writing xcmg and ycmg
c          write(sds_name,'(A4)') "xcmg"
c 	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
c          status=sfwdata(sds_id,start,stride,edges,txcmg)
c 	 write(6,*) "status sfwdata ",status
c          status = sfendacc(sds_id)
c          write(sds_name,'(A4)') "ycmg"
c 	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
c          status=sfwdata(sds_id,start,stride,edges,tycmg)
c 	 write(6,*) "status sfwdata ",status
c          status = sfendacc(sds_id)
c writing twv,twvi
	 
c        write(sds_name,'(A3)') "twv"
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,twv)
c	 write(6,*) "status sfwdata ",status
cc       status = sfendacc(sds_id)
         ihdf=0
         if (iout.ne.0) then
c writing aerosol algorithm bands (toa)

cc       write(sds_name,'(A7,A2)') "toaband",suffix(1)
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,aerob1)
cc       write(sds_name,'(A7,A2)') "toaband",suffix(2)
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,aerob2)
cc       write(sds_name,'(A7,A2)') "toaband",suffix(4)
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,aerob4)
cc       write(sds_name,'(A7,A2)') "toaband",suffix(5)
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,aerob5)
cc       write(sds_name,'(A7,A2)') "toaband",suffix(7)
c	 sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
cc       status=sfwdata(sds_id,start,stride,edges,aerob7)
	 
         write(sds_name,'(A4)') "twvi"
 	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,twvi)
 	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A4)') "tozi"
 	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tozi)
 	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A6)') "tpresi"
 	 sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,tp)
 	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
	 endif
         if  (iout.ne.0) then
	 write(sds_name,'(A11)') "taero-band1"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,taero3)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	 write(sds_name,'(A5)') "taero"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,taero)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	 write(sds_name,'(A5)') "tresi"
	sds_id=sfcreate(sd_id,sds_name,DFNT_FLOAT32,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tresi)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	   write(sds_name,'(A8)') "tratiob1"
	   write(6,*) "tratiob1 ",tratiob1(672,709),tratiob1(709,672)
	sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tratiob1)
	 write(6,*) "status sfwdata ",status
c	 status = sfendacc(sds_id)
	 scalefactor=1000.
	 offset=0.0
	 status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
	 status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
	 status = sfendacc(sds_id)
	   write(sds_name,'(A8)') "tratiob2"
	sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tratiob2)
	 write(6,*) "status sfwdata ",status
c	 status = sfendacc(sds_id)
	 scalefactor=1000.
	 offset=0.0
	 status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
	 status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
	 status = sfendacc(sds_id)
	 write(sds_name,'(A6)') "nbiter"
	sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tnit)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	   write(sds_name,'(A8)') "tratiob7"
	sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tratiob7)
	 write(6,*) "status sfwdata ",status
c	 status = sfendacc(sds_id)
	 scalefactor=1000.
	 offset=0.0
	 status=sfsattr(sds_id, "scale_factor", DFNT_FLOAT32, 1, scalefactor)
	 status=sfsattr(sds_id, "offset", DFNT_FLOAT32, 1, offset)
	 status = sfendacc(sds_id)
	 write(sds_name,'(A6)') "nbiter"
	sds_id=sfcreate(sd_id,sds_name,DFNT_INT16,rank,dim_sizes)
	 status=sfwdata(sds_id,start,stride,edges,tnit)
	write(6,*) "status sfwdata ",status
	 status = sfendacc(sds_id)
	 endif
	 
         write(sds_name,'(A5)') "CLOUD"
	 sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,cloud)
	 write(6,*) "status sfwdata ",status
	 write(6,*) "writing attribute data for CLOUD band ",sds_id
	 status=set_qamap(sds_id)
	 write(6,*) "status set qa attribute ",status
         status = sfendacc(sds_id)
         write(sds_name,'(A6)') "IPFLAG"
	 sds_id=sfcreate(sd_id,sds_name,DFNT_UINT8,rank,dim_sizes)
         status=sfwdata(sds_id,start,stride,edges,ipflag)
	 write(6,*) "status sfwdata ",status
         status = sfendacc(sds_id)
c writing projection metadata
	 xd=nc
	 yd=nr
	 uplx=x0+gsize*(frcl+col0-2)
	 uply=y0-gsize*(frln+row0-2)
	 lorx=uplx+nc*gsize
	 lory=uply-nr*gsize
	 zone=utmzone
	 sphere=12
	 
	 status=set_proj(sd_id,xd,yd,uplx,uply,lorx,lory,zone,sphere,wbc,ebc,nbc,sbc)
	 write(6,*) "status set projection attribute ",status
  	    
!closing the HDF file        	 
         status = sfend(sd_id) 	 
	 write(6,*) "status sfend ",status
	
	stop
	end
	

