cccccccccccccccccccc
c Modified on 11/16/2012 by Gail Schmidt, USGS ERO
c - Use snfindex to find the index of the specified SDS rather than relying
c   on hard-coded SDS index values
c
c Modified on 11/9/2012 by Gail Schmidt, USGS EROS
c - Use the fill_QA band to determine if the current pixel is fill
c - Don't set cloud-related QA bits to on if the current pixel is fill
c - Added error checking to file write/close tasks
c
c Modified on 10/16/2012 by Gail Schmidt, USGS EROS
c - pass in the pixel size from the file metadata vs. hard coding it to 28.5
c - modified to process individual cloud QA bands vs. the previous single
c   QA band with packed bits
cccccccccccccccccccc
         Program readl7sr
         
         character*200 filein
         integer ii
         integer sfstart, sfselect, sfrdata, sfendacc, sfend
         integer sfscatt, sfn2index, sfwdata, sfginfo
         integer sd_id, sds_id, sds_index, status
         integer start(5), edges(5), stride(5)
         integer DFACC_READ,DFACC_RDWR,DFNT_CHAR8
         parameter (DFACC_READ = 1)
         parameter (DFACC_RDWR = 3)
         integer nc,nr,mnc,mnr
         parameter (DFNT_CHAR8 = 4)
         parameter (mnc=9000)
         parameter (mnr=9000)
         integer(1) QA_ON, QA_OFF
         parameter (QA_OFF=INT(0,1))
         parameter (QA_ON=INT(255,1))
         integer ERROR, SUCCESS
         parameter (ERROR = 1)
         parameter (SUCCESS = 0)

         integer(1), allocatable :: cloud_qa(:,:)
         integer(1), allocatable :: cloud_shad_qa(:,:)
         integer(1), allocatable :: cloud_adja_qa(:,:)
         integer(1), allocatable :: snow_qa(:,:)
         integer(1), allocatable :: fill_qa(:,:)
         integer(1), allocatable :: tmpbit_qa(:,:)
         integer(2), allocatable :: temp(:,:)
         integer(2), allocatable :: band1(:,:)
         integer(2), allocatable :: band3(:,:)
         integer(2), allocatable :: band2(:,:)
         integer(2), allocatable :: band5(:,:)
         real tclear,ts,tv,fs,tcloud,cfac,cldh,pixsize,tna
         real fack,facl
         character*80 sds_name
         integer rank,data_type
         integer n_attrs
         integer dim_sizes(5)
         integer mband5,mband5k,mband5l
         double precision mclear
         real pclear
         integer*8 nbclear,nbval,nbcloud
         real t6
         integer cldhmin,cldhmax
         
         cfac=6.
         dtr=atan(1.)/45.
         pixsize=28.5 
         
         call getarg(1,filein)
         ii=index(filein," ")-1
         write(6,*) "temp clear [k] , ts,tv,fs,truenorthadj,pixsize"
         read(5,*)  tclear,ts,tv,fs,tna,pixsize
         write(6,*)  tclear,ts,tv,fs,tna,pixsize
         tclear=tclear-273.0
                 
cc Read from the HDF file 
         sd_id= sfstart(filein(1:ii),DFACC_RDWR)
         
         write(6,*) "sd_id", sd_id
         if (sd_id.lt.0) then
         write(6,'(A25,1X,A80)')
     &  "WW cannot open input ",filein
         call exit(ERROR)
         endif

c reading qa       
c allocate memory for qa
       sds_index = sfn2index(sd_id, "cloud_QA")
       if (sds_index .lt. 0) then
           write(6,*) "error finding cloud qa data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "cloud sds_id", sds_id
       status= sfginfo(sds_id, sds_name, rank, dim_sizes, data_type,
     s  n_attrs)
       if (status .ne. 0) then
           write(6,*) "error finding cloud qa data"
           call exit(ERROR)
       endif
       write(6,*) "sdsname ",sds_name
       write(6,*) dim_sizes(1),dim_sizes(2)
       nc= dim_sizes(1)
       nr=dim_sizes(2)
       allocate (cloud_qa(nc,nr),stat=ierr)
       allocate (cloud_shad_qa(nc,nr),stat=ierr)
       allocate (cloud_adja_qa(nc,nr),stat=ierr)
       allocate (snow_qa(nc,nr),stat=ierr)
       allocate (fill_qa(nc,nr),stat=ierr)
       allocate (tmpbit_qa(nc,nr),stat=ierr)
       allocate (temp(nc,nr),stat=ierr)
       allocate (band1(nc,nr),stat=ierr)
       allocate (band2(nc,nr),stat=ierr)
       allocate (band3(nc,nr),stat=ierr)
       allocate (band5(nc,nr),stat=ierr)

cc read cloud qa data 
       start(1)=0
       start(2) = 0
       edges(1) = nc
       edges(2) = nr
       stride(1) = 1
       stride(2) = 1
       status = sfrdata(sds_id, start, stride, edges, cloud_qa)
       if (status .ne. 0) then
           write(6,*) "error reading cloud qa data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading cloud shadow qa data
       sds_index = sfn2index(sd_id, "cloud_shadow_QA")
       if (sds_index .lt. 0) then
           write(6,*) "error finding shadow qa data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       if (status .ne. 0) then
           write(6,*) "error finding shadow qa data"
           call exit(ERROR)
       endif
c       write(6,*) "cloud shadow sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, cloud_shad_qa)
       if (status .ne. 0) then
           write(6,*) "error reading shadow qa data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading adjacent cloud qa data
       sds_index = sfn2index(sd_id, "adjacent_cloud_QA")
       if (sds_index .lt. 0) then
           write(6,*) "error finding adj cloud qa data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       if (status .ne. 0) then
           write(6,*) "error finding adj cloud qa data"
           call exit(ERROR)
       endif
c       write(6,*) "adjacent cloud sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, cloud_adja_qa)
       if (status .ne. 0) then
           write(6,*) "error reading adj cloud qa data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading snow qa data
       sds_index = sfn2index(sd_id, "snow_QA")
       if (sds_index .lt. 0) then
           write(6,*) "error finding snow qa data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       if (status .ne. 0) then
           write(6,*) "error finding snow qa data"
           call exit(ERROR)
       endif
c       write(6,*) "snow qa sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, snow_qa)
       if (status .ne. 0) then
           write(6,*) "error reading snow qa data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading fill qa data
       sds_index = sfn2index(sd_id, "fill_QA")
       if (sds_index .lt. 0) then
           write(6,*) "error finding fill qa data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       if (status .ne. 0) then
           write(6,*) "error finding fill qa data"
           call exit(ERROR)
       endif
c       write(6,*) "fill qa sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, fill_qa)
       if (status .ne. 0) then
           write(6,*) "error reading fill qa data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading for tmpbit is not necessary since it doesn't exist in the lndsr
c product and it is set to zero below before using

c reading temperature       
       sds_index = sfn2index(sd_id, "band6")
       if (sds_index .lt. 0) then
           write(6,*) "error reading temperature data"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "temperature sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, temp)
       if (status .ne. 0) then
           write(6,*) "error reading temperature data"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading band1 (blue)
       sds_index = sfn2index(sd_id, "band1")
       if (sds_index .lt. 0) then
           write(6,*) "error reading blue band"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "band 1 sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, band1)
       if (status .ne. 0) then
           write(6,*) "error reading blue band"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading band2 (green)
       sds_index = sfn2index(sd_id, "band2")
       if (sds_index .lt. 0) then
           write(6,*) "error reading green band"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "band 2 sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, band2)
       if (status .ne. 0) then
           write(6,*) "error reading green band"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading band3 (red)
       sds_index = sfn2index(sd_id, "band3")
       if (sds_index .lt. 0) then
           write(6,*) "error reading red band"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "band 3 sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, band3)
       if (status .ne. 0) then
           write(6,*) "error reading red band"
           call exit(ERROR)
       endif
c       write(6,*) "status", status
c reading band5 (swir)
       sds_index = sfn2index(sd_id, "band5")
       if (sds_index .lt. 0) then
           write(6,*) "error reading swir band"
           call exit(ERROR)
       endif
       sds_id    = sfselect(sd_id, sds_index)
       write(6,*) "band 5 sds_id", sds_id
       status = sfrdata(sds_id, start, stride, edges, band5)
       if (status .ne. 0) then
           write(6,*) "error reading swir band"
           call exit(ERROR)
       endif
c       write(6,*) "status", status

c printing some values for checking out
c       i=3953
c       j=3153
c       write(6,*) " i,j ", i,j
c       write(6,*) "cloud_qa(j,i) ",cloud_qa(j,i)
c       write(6,*) "cloud_shad_qa(j,i) ",cloud_shad_qa(j,i)
c       write(6,*) "cloud_adja_qa(j,i) ",cloud_adja_qa(j,i)
c       write(6,*) "snow_qa(j,i) ",snow_qa(j,i)
c       write(6,*) "fill_qa(j,i) ",fill_qa(j,i)
c       write(6,*) "{garbage} tmpbit_qa(j,i) ",tmpbit_qa(j,i)
c       write(6,*) "temp(j,i) ",temp(j,i)
c       write(6,*) "band1(j,i) ",band1(j,i)
c       write(6,*) "band3(j,i) ",band3(j,i)
c       write(6,*) "band5(j,i) ",band5(j,i)
       
c note i is the line and j is the column  
       write(6,*) "resetting cloud,cloudadja,cloudshad bits,tmpbit"
       do i=1,nr
       do j=1,nc
c reset cloud,cloud shadow,adjacent cloud bit
       cloud_qa(j,i) = QA_OFF
       cloud_adja_qa(j,i) = QA_OFF
       cloud_shad_qa(j,i) = QA_OFF
       tmpbit_qa(j,i) = QA_OFF
       enddo
       enddo
       
c update the cloud mask
       write(6,*) "updating cloud mask"
       nbclear=0
       mclear=0.
       nbval=0
       nbcloud=0
      
c compute the average temperature of the clear data       
       do i=1,nr
       do j=1,nc
       if (fill_qa(j,i).eq.QA_OFF) then
       nbval=nbval+1
       anom=band1(j,i)-band3(j,i)/2.
       t6=temp(j,i)/100.
       if (snow_qa(j,i).eq.QA_ON) then
       continue
       else
       if ((anom.gt.300.).and.(band5(j,i).gt.300)
     s   .and.(t6.lt.tclear)) then
       cloud_qa(j,i)=QA_ON
       nbcloud=nbcloud+1
       else
       if ((band1(j,i).gt.3000).and.(t6.lt.tclear)) then
       cloud_qa(j,i)=QA_ON
       nbcloud=nbcloud+1
       else
       mclear=mclear+t6/10000.
       nbclear=nbclear+1
       endif
       endif
       endif
       endif
       enddo
       enddo
c 
        if (nbclear.gt.0) then
        mclear=mclear*10000./nbclear
        endif
        pclear=(nbclear*100./nbval)
        write(6,*) nbclear-nbval,nbcloud
        write(6,*) "average clear temperature  %clear", mclear,pclear
        if (pclear.gt.5.) then
        tclear=mclear       
        endif
             
c update the adjacent cloud bit. only set for non-fill values
       write(6,*) "updating adjacent cloud bit"
       do i=1,nr
       do j=1,nc
       if (cloud_qa(j,i).eq.QA_ON) then
c       write(6,*) "i found a cloud at ",i,j,"with t=",temp(j,i)/100.
       do k=i-5,i+5
       do l=j-5,j+5
       if ((k.ge.1).and.(k.le.nr).and.(l.ge.1).and.(l.le.nc)) then
       if ((cloud_adja_qa(l,k).eq.QA_ON).or.(cloud_qa(l,k).eq.QA_ON)
     &     .or.(fill_qa(l,k) .eq. QA_ON)) then
       continue
       else
       cloud_adja_qa(l,k)=QA_ON
       endif
       endif
       enddo
       enddo
       endif
       enddo
       enddo
       
C compute cloud shadow  
       write(6,*) "looking for cloud shadow"
       fack=cos(fs*dtr)*tan(ts*dtr)/pixsize
       facl=sin(fs*dtr)*tan(ts*dtr)/pixsize
       do i=1,nr
       do j=1,nc
c       if ((i.eq.5948).and.(j.eq.5291)) then
c       write(6,*) "we are at 5948,5291"
c       endif
       if (cloud_qa(j,i).eq.QA_ON) then
           tcloud=temp(j,i)/100.
           cldh=(tclear-tcloud)*1000./cfac
           if (cldh.lt.0.) cldh=0.
           cldhmin=cldh-1000.
           cldhmax=cldh+1000.
           mband5=9999
           do icldh=cldhmin/10,cldhmax/10
            cldh=icldh*10.
            k=i+fack*cldh
            l=j-facl*cldh
            if ((k.ge.1).and.(k.le.nr).and.(l.ge.1)
     &            .and.(l.le.nc)) then
            if ((band5(l,k).lt.800).and.
     &            ((band2(l,k)-band3(l,k)).lt.100)) then
            if ((cloud_adja_qa(l,k).eq.QA_ON).or.
     &            (cloud_qa(l,k).eq.QA_ON).or.
     &            (cloud_shad_qa(l,k).eq.QA_ON).or.
     &            (fill_qa(l,k).eq.QA_ON)) then
             continue
             else
c store the value of band5 as well as the l and k value
             if (band5(l,k).lt.mband5) then
             mband5=band5(l,k)
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
             cloud_shad_qa(l,k)=QA_ON
             endif
         endif
         enddo
         enddo
c dilate the cloud shadowc
       write(6,*) "dilating cloud shadow"
       do i=1,nr
       do j=1,nc
       if (cloud_shad_qa(j,i).eq.QA_ON) then
       do k=i-3,i+3
       do l=j-3,j+3
       if ((k.ge.1).and.(k.le.nr).and.(l.ge.1).and.(l.le.nc)) then
       if ((cloud_adja_qa(l,k).eq.QA_ON).or.(cloud_qa(l,k).eq.QA_ON)
     &   .or.(cloud_shad_qa(l,k).eq.QA_ON).or.(fill_qa(l,k).eq.QA_ON))
     &   then
       continue
       else
       if (tmpbit_qa(l,k).eq.QA_ON) then
       continue
       else
       tmpbit_qa(l,k)=QA_ON
       endif
       endif
       endif
       enddo
       enddo
       endif
       enddo
       enddo
c update the cloud shadow and clear QA for fill pixels
       write(6,*) "updating cloud shadow"
       do i=1,nr
       do j=1,nc
       if (tmpbit_qa(j,i).eq.QA_ON) then
       cloud_shad_qa(j,i)=QA_ON
       tmpbit_qa(j,i)=QA_OFF
       endif

c      clear the QA info for the fill pixels
       if (fill_qa(j,i) .eq. QA_ON) then
           cloud_qa(j,i) = QA_OFF
           cloud_adja_qa(j,i) = QA_OFF
           cloud_shad_qa(j,i) = QA_OFF
       endif 
       enddo
       enddo
       
c updating the cloud qa and closing SDS
        sds_index = sfn2index(sd_id, "cloud_QA")
        if (sds_index .lt. 0) then
            write(6,*) "error finding cloud qa data"
            call exit(ERROR)
        endif
        sds_id= sfselect(sd_id,sds_index)
        status=sfwdata(sds_id, start, stride,edges,cloud_qa)
        if (status .ne. 0) then
            write(6,*) "error writing cloud qa. status:", status
            call exit(ERROR)
        endif
c        write(6,*) "write cloud qa status", status
        status = sfendacc(sds_id)
        if (status .ne. 0) then
            write(6,*) "error closing cloud qa. status:", status
            call exit(ERROR)
        endif
c        write(6,*) "close cloud qa status", status
c updating the cloud shadow qa and closing SDS
        sds_index = sfn2index(sd_id, "cloud_shadow_QA")
        if (sds_index .lt. 0) then
            write(6,*) "error writing shadow qa. status:", status
            call exit(ERROR)
        endif
        sds_id= sfselect(sd_id,sds_index)
        status=sfwdata(sds_id, start, stride,edges,cloud_shad_qa)
        if (status .ne. 0) then
            write(6,*) "error writing shadow qa. status:", status
            call exit(ERROR)
        endif
c        write(6,*) "write cloud shadow qa status", status
        status = sfendacc(sds_id)
        if (status .ne. 0) then
            write(6,*) "error closing shadow qa. status:", status
            call exit(ERROR)
        endif
c        write(6,*) "close cloud shadow qa status", status
c updating the adjacent cloud qa and closing SDS
        sds_index = sfn2index(sd_id, "adjacent_cloud_QA")
        if (sds_index .lt. 0) then
            write(6,*) "error writing adjacent cloud qa. status:",status
            call exit(ERROR)
        endif
        sds_id= sfselect(sd_id,sds_index)
        status=sfwdata(sds_id, start, stride,edges,cloud_adja_qa)
        if (status .ne. 0) then
            write(6,*) "error writing adjacent cloud qa. status:",status
            call exit(ERROR)
        endif
c        write(6,*) "write adjacent cloud qa status", status
        status = sfendacc(sds_id)
        if (status .ne. 0) then
            write(6,*) "error closing adjacent cloud qa. status:",status
            call exit(ERROR)
        endif
c        write(6,*) "close adjacent cloud qa status", status
c update the global attribute to reflect cloud processing
        status=sfscatt(sd_id, "Cloud Mask Algo Version ",DFNT_CHAR8,
     s  22, "CMReflectanceBasedv1.0")
        if (status .ne. 0) then
            write(6,*) "error writing global attributes"
            call exit(ERROR)
        endif
c        write(6,*) "status write attribute", status
c close HDF file
        status = sfend(sd_id)
        if (status .ne. 0) then
            write(6,*) "error closing surface reflectance file"
            call exit(ERROR)
        endif
        deallocate(cloud_qa,cloud_shad_qa,cloud_adja_qa,snow_qa,
     &   tmpbit_qa)
       deallocate(temp,band1,band2,band3,band5)
       call exit(SUCCESS)
       end
