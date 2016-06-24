        subroutine subaeroret(iband1,iband3,xts,xtv,xfi,pres,uoz,uwv,erelc,troatm,
     c       tpres,aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,iv,raot,residual,nit,snext)
        IMPLICIT NONE
c driver for the atmopsheric correction program
c read up the look up table and perform correction
c look up table variable
c   Aerosol model (string) and length of the string
c local variable
        real rotoa,xts,xtv,xfi,raot550nm,uoz,uwv,pres
	real residual
	integer ib,iv
	real roslamb,ros1,ros3,raot1,raot2
        integer band,flagn
       
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
       real aot550nm(22),xtest
       real xaot,next,snext
       real tpres(7)
       integer retval
       integer iaot
        real tgo,roatm,ttatmg,satm,xrorayp
	real troatm(16),roref,trosur(16),thinv(16)
	real aratio1,pratio,aratio2,eratio,inpaot,eaot,th1,th3
	real peratio,pros1,pros3
        real raot
        real erelc(8)
	integer iband1,iband3,iband,nit,iter
	character*80 err_msg
	integer nb
	
	   
	
c correct band 3 and band 1 with increasing AOT (using pre till ratio is equal to erelc(3)
        iaot=1
	pratio=erelc(iband3)/erelc(iband1)
c pratio is the targeted ratio between the surface reflectance in band 4 (ros4) and 2 (ros2)
	aratio1=1000.
	aratio2=2000.
	ros1=1.0
	ros3=1.0
	raot1=0.0001
	flagn=0
	th1=1.E-02
	th3=1.E-02
	nit=0
        do while ((iaot.le.22).and.(aratio1.gt.pratio).and.(ros1.gt.th1).and.(ros3.gt.th3).and.((aratio1-0.01).lt.aratio2)
     s  	.and.(nit.lt.30))
c the ratio decrease as the AOT increase the exit condition in this loop are when two values of AOT can be found that bracket the predicted ratio pratio
        nit=nit+1	
	ros1=-1.
	ros3=-1.
c if flagn is set... we start converge to the AOT bounds by dichotomy to increase the accuracy of the retrieval	
	if (flagn.eq.0) then
	raot550nm=aot550nm(iaot)
	else
	raot550nm=(raot1+aot550nm(iaot))/2.
	endif
	
	iter=0
	do while ((ros1.lt.th1).or.(ros3.lt.th3)) 	
	if (iter.gt.0) then 
	if (iaot.ge.2) then
	raot550nm=(raot550nm+aot550nm(iaot-1))/2.
	else
	if (iv.eq.1) write(6,*) "report inversion failed "
	goto 91
	endif
	endif
	iband=iband3
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros3=roslamb
	iband=iband1
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros1=roslamb
	iter=iter+1
	if (iv.eq.1) write(6,*) "report ",nit, raot550nm,ros1,ros3,aratio1,pratio
	enddo
	
	if ((iter.gt.1).or.(flagn.eq.1)) then
	flagn=1
	else
	iaot=iaot+1
	endif
	
	if ((ros1.gt.th1).and.(ros3.gt.th3)) then
	aratio2=aratio1
	raot2=raot1
	raot1=raot550nm
	aratio1=ros3/ros1
	if (iv.eq.1) write(6,*) "report ",nit, raot550nm,ros1,ros3,aratio1,pratio
	else
	if (iv.eq.1) write(6,*) "report negative",raot550nm,ros1,ros3
        endif 
	enddo
c once the two value of AOT  (raot2 and raot1) that gives ratios that bracket the predicted ratio are found
C they are used to estimate the aot (eaot) using linear interpolation	
        if ((aratio1.gt.pratio).and.(aratio2.gt.pratio)) then
        if (raot1.lt.raot2) then 
	raot550nm=raot1
	else
	raot550nm=raot2
	endif
	goto 91
	endif	
	
	eaot=(aratio1-pratio)*(raot2-raot1)/(aratio1-aratio2)+raot1

c this eaot is refined by recomputing by performing an additional iteration
	
	raot550nm=eaot
	iband=iband3
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros3=roslamb
	iband=iband1
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros1=roslamb
	
c	write(6,*) "report ",ros1,ros3,eaot,ros3/ros1,pratio,raot1,raot2
	eratio=ros3/ros1
	if (abs(eratio-aratio1).gt.abs(eratio-aratio2)) then
	raot2=eaot
	aratio2=eratio
	eaot=(aratio1-pratio)*(raot2-raot1)/(aratio1-aratio2)+raot1
	else
	raot1=eaot
	aratio1=eratio
	eaot=(aratio1-pratio)*(raot2-raot1)/(aratio1-aratio2)+raot1
	endif
	
	raot550nm=eaot
	iband=iband3
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros3=roslamb
	iband=iband1
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros1=roslamb
c the last step of the algorithm is done by making very small increase or decrease of the estimated aot (0.0005) 
c the value of eaot giving the closest ratio to pratio is finally selected

        eratio=ros3/ros1
	raot550nm=eaot
	if (raot550nm.lt.0.01) goto 91
	peratio=1000.
	if (eratio.gt.pratio) then
	do while ((eratio.gt.pratio).and.(peratio.gt.eratio))
	pros1=ros1
	pros3=ros3
	raot550nm=raot550nm+0.005
c	write(6,*) "I am here ",eratio,pratio
	iband=iband3
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros3=roslamb
	iband=iband1
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros1=roslamb
	peratio=eratio
	eratio=ros3/ros1
c	write(6,*) eratio,pratio,ros1,ros3
	enddo
	if (abs(eratio-pratio).gt.abs(peratio-pratio)) then
	    raot550nm=raot550nm-0.005
	    eratio=peratio
	    ros1=pros1
	    ros3=pros3
	    endif


	else
	
	peratio=0.
	do while ((eratio.lt.pratio).and.(peratio.lt.eratio))
	pros1=ros1
	pros3=ros3
	raot550nm=raot550nm-0.005
	iband=iband3
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros3=roslamb
	iband=iband1
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
        ros1=roslamb
	peratio=eratio
	eratio=ros3/ros1
	enddo
	if (abs(eratio-pratio).gt.abs(peratio-pratio)) then
	    raot550nm=raot550nm+0.005
	    eratio=peratio
	    ros1=pros1
	    ros3=pros3
	    endif
	
        endif
c compute model residual
  91   residual=abs(ros3-ros1*pratio)
        nb=1
c       write(6,*) "INITIAL RESIDUAL ",residual
       do iband=1,6
       if (erelc(iband).gt.0.) then
	call atmcorlamb2(xts,xtv,xfi,raot550nm,iband,pres,tpres,
     s       aot550nm,rolutt,
     s       transt,xtsstep,xtsmin,xtvstep,xtvmin,
     s       sphalbt,normext,
     s       tsmax,tsmin,nbfic,nbfi,tts,indts,ttv,
     s       uoz,uwv,tauray,
     s       ogtransa0,ogtransa1,ogtransb0,ogtransb1,
     s       ogtransc0,ogtransc1,
     s       wvtransa,wvtransb,wvtransc,oztransa,     
     s       troatm(iband),roslamb,tgo,roatm,ttatmg,satm,xrorayp,next,
     s       err_msg,retval)
c       	ros=roslamb
        if (iband.eq.iband3) then
	 snext=next
	else
	
	residual=residual+abs(roslamb-ros1*(erelc(iband)/erelc(iband1)))
	nb=nb+1
         endif
c	 write(6,*) "band ",iband, "residual ", abs(roslamb-ros1*(erelc(iband)/erelc(iband1)))
	endif
	enddo
	residual=residual/(nb-1)
c	write(6,*) "MODEL RESIDUAL ",residual
	raot=raot550nm
	
	
 99	continue
c         write(6,*) "finalreport ",raot550nm,eratio,pratio,troatm(iband3),troatm(iband1),ros3,ros1,"residual ",residual
	
	raot=raot550nm
  	return
	end
	
	
	
