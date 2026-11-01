! Scenario adapter for the pinned SCALE build. Native collision/drag/growth laws
! are retained. This file is inserted into scale_atmos_phy_mp_sdm by the builder.
! STATE
  logical, save :: sl_ready=.false.
  integer, save :: sl_mode=0, sl_layers, sl_packets, sl_calls=0, sl_au, sl_lu
  integer(kind=8), save :: sl_rng
  real(RP), save :: sl_depth, sl_dx, sl_flux, sl_rate, sl_end, sl_gather
  real(RP), save :: sl_c,sl_b,sl_dmin,sl_dmax,sl_time=0,sl_h
  real(RP), save :: sl_in=0,sl_dep=0,sl_supply=0,sl_out=0,sl_vinitial=0
  real(RP), save :: sl_maxerr=0,sl_maxd=0,sl_peak=0,sl_rounding=0,sl_removed_heat=0
  real(RP), allocatable, save :: sl_temp(:),sl_dry(:),sl_vap(:),sl_carry(:,:)
  real(RP), allocatable, save :: sl_oldre(:),sl_oldrp(:),sl_oldrho(:)
! PROCEDURES
  function sl_random() result(x)
    real(RP) :: x
    sl_rng=modulo(16807_8*sl_rng,2147483647_8)
    x=real(sl_rng,RP)/2147483647.0_RP
  end function

  subroutine sl_setgas(d,q,rhot,temp,vap,dry)
    use scale_const, only: rd=>CONST_Rdry,rv=>CONST_Rvap,cp=>CONST_CPdry,p0=>CONST_PRE00
    use scale_atmos_thermodyn, only: CPw=>AQ_CP
    real(RP),intent(out)::d,q,rhot
    real(RP),intent(in)::temp,vap,dry
    real(RP)::r,cpair,p
    d=dry+vap
    q=vap/d
    r=rd*(1-q)+rv*q
    cpair=cp*(1-q)+CPw(I_QV)*q
    p=d*r*temp
    rhot=d*temp*(p0/p)**(r/cpair)
  end subroutine

  subroutine sl_start(dens,rhot,qtrc,sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_rk,sd_r,sd_asl,sdi,active)
    use scale_const, only: rd=>CONST_Rdry,rv=>CONST_Rvap,cp=>CONST_CPdry
    use scale_tracer, only: QA
    use scale_atmos_thermodyn, only: CPw=>AQ_CP
    use m_sdm_fluidconv, only: sdm_rhot_qtrc2p_t
    real(RP),intent(inout)::dens(KA,IA,JA),rhot(KA,IA,JA),qtrc(KA,IA,JA,QA)
    integer,intent(in)::sd_num
    integer(DP),intent(inout)::sd_n(sd_num)
    integer(i2),intent(inout)::sd_liqice(sd_num)
    real(RP),intent(inout)::sd_x(sd_num),sd_y(sd_num),sd_rk(sd_num),sd_r(sd_num),sd_asl(:,:)
    type(sdicedef),intent(inout)::sdi
    logical,intent(out)::active
    integer::fu,k,i,j,l,n,p,slot,kk
    integer(DP)::countn,nn,leftn
    real(RP)::tt,pp,ee,vap,dry,add,es,tlocal,d,m,pi,requested
    real(RP)::pres(KA,IA,JA),temp(KA,IA,JA)
    logical::exists
    pi=acos(-1.0_RP)
    if(.not.sl_ready)then
      sl_ready=.true.
      inquire(file='lunar_forcing.txt',exist=exists)
      if(exists)then
        open(newunit=fu,file='lunar_forcing.txt',status='old')
        read(fu,*)sl_mode,sl_layers,sl_packets,sl_rng
        read(fu,*)sl_depth,sl_dx,sl_flux,sl_rate,sl_end,sl_gather
        read(fu,*)sl_c,sl_b,sl_dmin,sl_dmax
        if(sl_layers/=KE-KS .or. sl_mode<1 .or. sl_mode>3) error stop 'Invalid lunar forcing grid/mode'
        allocate(sl_temp(KA),sl_dry(KA),sl_vap(KA),sl_carry(IA,JA))
        allocate(sl_oldre(sd_num),sl_oldrp(sd_num),sl_oldrho(sd_num))
        sl_carry=0
        do l=1,sl_layers
          k=KS+l-1
          read(fu,*)tt,pp,ee
          sl_temp(k)=tt
          sl_vap(k)=ee/(rv*tt)
          sl_dry(k)=(pp-ee)/(rd*tt)
        enddo
        close(fu)
        do k=1,KA
          kk=max(KS,min(KS+sl_layers-1,k))
          sl_temp(k)=sl_temp(kk);sl_dry(k)=sl_dry(kk);sl_vap(k)=sl_vap(kk)
          do j=1,JA
          do i=1,IA
            call sl_setgas(dens(k,i,j),qtrc(k,i,j,I_QV),rhot(k,i,j),sl_temp(k),sl_vap(k),sl_dry(k))
          enddo
          enddo
        enddo
        sl_vinitial=sum(sl_vap(KS:KS+sl_layers-1))*sl_depth/sl_layers
        open(newunit=sl_au,file='lunar_arrivals.bin',status='new',access='stream',form='unformatted')
        open(newunit=sl_lu,file='lunar_ledger.txt',status='new')
        write(sl_lu,'(A)') '# time injected deposited vapor_source precipitation suspended vapor residual max_abs_residual Tmin Tmax maximum_D active peak_active removed_heat_J_m2'
      endif
    endif
    active=.true.
    if(sl_mode==0)return
    sl_h=min(sdm_dtcmph(3),sl_end-sl_time)
    if(sl_h<1e-8_RP)then
      active=.false.;return
    endif
    call sdm_rhot_qtrc2p_t(rhot,qtrc,dens,pres,temp)
    do k=KS,KS+sl_layers-1
    do j=JS,JE
    do i=IS,IE
      if(sl_mode==1)then
        call sl_setgas(dens(k,i,j),qtrc(k,i,j,I_QV),rhot(k,i,j),sl_temp(k),sl_vap(k),sl_dry(k))
      else
        vap=dens(k,i,j)*qtrc(k,i,j,I_QV)
        dry=dens(k,i,j)-vap
        tlocal=temp(k,i,j)
        if(sl_mode==3)then
          ! Diagnostic sensible cooling needed to maintain the imposed profile.
          sl_removed_heat=sl_removed_heat+(dry*cp+vap*CPw(I_QV))*(tlocal-sl_temp(k))*sl_depth/sl_layers/real((IE-IS+1)*(JE-JS+1),RP)
          tlocal=sl_temp(k)
        endif
        ! Same Murphy--Koop liquid saturation curve as the comparison column.
        es=exp(54.842763_RP-6763.22_RP/tlocal-4.210_RP*log(tlocal)+0.000367_RP*tlocal &
          +tanh(0.0415_RP*(tlocal-218.8_RP))*(53.878_RP-1331.22_RP/tlocal-9.44523_RP*log(tlocal)+0.014025_RP*tlocal))
        add=min(sl_rate*sl_h,max(0.0_RP,es/(rv*tlocal)-vap))
        sl_supply=sl_supply+add*sl_depth/sl_layers/real((IE-IS+1)*(JE-JS+1),RP)
        call sl_setgas(dens(k,i,j),qtrc(k,i,j,I_QV),rhot(k,i,j),tlocal,vap+add,dry)
      endif
    enddo
    enddo
    enddo
    ! Independent seed stream; changing collision random draws does not change it.
    slot=1
    do j=JS,JE
    do i=IS,IE
      requested=sl_flux*sl_dx**2*sl_h+sl_carry(i,j)
      countn=int(requested,kind=DP)
      sl_carry(i,j)=requested-real(countn,RP)
      leftn=countn
      do p=1,sl_packets
        nn=leftn/int(sl_packets-p+1,DP)
        leftn=leftn-nn
        if(nn==0)cycle
        do while(slot<=sd_num)
          if(sd_rk(slot)<PREC2INVALID)exit
          slot=slot+1
        enddo
        if(slot>sd_num)error stop 'Lunar seed allocation exhausted; increase capacity'
        n=slot;slot=slot+1
        d=exp(log(sl_dmin)+sl_random()*log(sl_dmax/sl_dmin))
        m=sl_c*d**sl_b
        sd_n(n)=nn;sd_liqice(n)=STAT_ICE
        sd_x(n)=(real(i-IS,RP)+0.5_RP)*sl_dx
        sd_y(n)=(real(j-JS,RP)+0.5_RP)*sl_dx
        sd_rk(n)=real(KS-1+sl_layers,RP)-1e-7_RP
        sd_r(n)=0;sd_asl(n,:)=0
        sdi%re(n)=d/2;sdi%rp(n)=d/2
        sdi%rho(n)=m/(4.0_RP*pi/3*(d/2)**3)
        sdi%tf(n)=230;sdi%mrime(n)=0;sdi%nmono(n)=1
        sl_in=sl_in+real(nn,RP)*m/(sl_dx**2*real((IE-IS+1)*(JE-JS+1),RP))
      enddo
    enddo
    enddo
  end subroutine

  subroutine sl_before_growth(sdi)
    type(sdicedef),intent(in)::sdi
    if(sl_mode==0)return
    sl_oldre=sdi%re;sl_oldrp=sdi%rp;sl_oldrho=sdi%rho
  end subroutine

  subroutine sl_limit_growth(dens,qtrc,sd_num,sd_n,sd_liqice,sd_ri,sd_rj,sd_rk,sdi,h)
    use scale_tracer, only: QA
    real(RP),intent(in)::dens(KA,IA,JA),qtrc(KA,IA,JA,QA),h
    integer,intent(in)::sd_num
    integer(DP),intent(in)::sd_n(sd_num)
    integer(i2),intent(in)::sd_liqice(sd_num)
    real(RP),intent(in)::sd_ri(sd_num),sd_rj(sd_num),sd_rk(sd_num)
    type(sdicedef),intent(inout)::sdi
    real(RP)::positive(KA,IA,JA),negative(KA,IA,JA),alpha(KA,IA,JA)
    real(RP)::m0,m1,dm,fac,desired,trial,vol,pi,area,available
    integer::n,k,i,j
    if(sl_mode==0)return
    pi=acos(-1.0_RP);vol=sl_dx**2*sl_depth/sl_layers
    area=sl_dx**2*real((IE-IS+1)*(JE-JS+1),RP)
    positive=0;negative=0;alpha=1
    do n=1,sd_num
      if(sd_rk(n)<VALID2INVALID .or. sd_liqice(n)/=STAT_ICE)cycle
      k=floor(sd_rk(n))+1;i=floor(sd_ri(n))+1;j=floor(sd_rj(n))+1
      m0=4*pi/3*sl_oldre(n)**2*sl_oldrp(n)*sl_oldrho(n)
      m1=4*pi/3*sdi%re(n)**2*sdi%rp(n)*sdi%rho(n)
      dm=real(sd_n(n),RP)*(m1-m0)
      positive(k,i,j)=positive(k,i,j)+max(dm,0.0_RP)
      negative(k,i,j)=negative(k,i,j)+min(dm,0.0_RP)
    enddo
    do k=KS,KS+sl_layers-1
    do j=JS,JE
    do i=IS,IE
      if(sl_mode==1)then
        available=sl_rate*h*vol
      else
        available=max(0.0_RP,dens(k,i,j)*qtrc(k,i,j,I_QV)*vol*(1-1e-12_RP)-negative(k,i,j))
      endif
      if(positive(k,i,j)>available)alpha(k,i,j)=available/positive(k,i,j)
    enddo
    enddo
    enddo
    do n=1,sd_num
      if(sd_rk(n)<VALID2INVALID .or. sd_liqice(n)/=STAT_ICE)cycle
      k=floor(sd_rk(n))+1;i=floor(sd_ri(n))+1;j=floor(sd_rj(n))+1
      m0=4*pi/3*sl_oldre(n)**2*sl_oldrp(n)*sl_oldrho(n)
      m1=4*pi/3*sdi%re(n)**2*sdi%rp(n)*sdi%rho(n)
      fac=alpha(k,i,j)
      if(m1>m0 .and. fac<1)then
        desired=m0+fac*(m1-m0)
        ! Interpolate the native morphology increment, then enforce its exact
        ! admitted mass with a common axis scale. This is a forcing adapter;
        ! timestep sensitivity is checked independently of SDM sampling.
        sdi%re(n)=sl_oldre(n)+fac*(sdi%re(n)-sl_oldre(n))
        sdi%rp(n)=sl_oldrp(n)+fac*(sdi%rp(n)-sl_oldrp(n))
        sdi%rho(n)=sl_oldrho(n)+fac*(sdi%rho(n)-sl_oldrho(n))
        trial=4*pi/3*sdi%re(n)**2*sdi%rp(n)*sdi%rho(n)
        fac=(desired/trial)**(1.0_RP/3)
        sdi%re(n)=sdi%re(n)*fac;sdi%rp(n)=sdi%rp(n)*fac
        m1=4*pi/3*sdi%re(n)**2*sdi%rp(n)*sdi%rho(n)
      endif
      sl_dep=sl_dep+real(sd_n(n),RP)*(m1-m0)/area
    enddo
  end subroutine

  subroutine sl_finish(dens,rhot,qtrc,sd_num,sd_n,sd_liqice,sd_r,sd_rk,sd_x,sd_y,sdi)
    use scale_tracer, only: QA
    use m_sdm_fluidconv,only:sdm_rhot_qtrc2p_t
    real(RP),intent(in)::dens(KA,IA,JA),rhot(KA,IA,JA),qtrc(KA,IA,JA,QA)
    integer,intent(in)::sd_num
    integer(DP),intent(in)::sd_n(sd_num)
    integer(i2),intent(in)::sd_liqice(sd_num)
    real(RP),intent(in)::sd_r(sd_num),sd_rk(sd_num),sd_x(sd_num),sd_y(sd_num)
    type(sdicedef),intent(in)::sdi
    real(RP)::suspended,vapor,area,m,d,weight,res,pi,rho,aspect,cell,maxd,source
    real(RP)::pres(KA,IA,JA),temp(KA,IA,JA)
    integer::n,live
    if(sl_mode==0)return
    pi=acos(-1.0_RP);area=sl_dx**2*real((IE-IS+1)*(JE-JS+1),RP)
    suspended=0;live=0;maxd=0
    sl_time=sl_time+sl_h;sl_calls=sl_calls+1
    do n=1,sd_num
      if(sd_rk(n)<PREC2INVALID)cycle
      if(sd_liqice(n)==STAT_ICE)then
        m=4*pi/3*sdi%re(n)**2*sdi%rp(n)*sdi%rho(n)
        d=2*max(sdi%re(n),sdi%rp(n));rho=sdi%rho(n);aspect=sdi%rp(n)/sdi%re(n)
      else
        m=4*pi/3*sd_r(n)**3*1000.0_RP
        d=2*sd_r(n);rho=1000;aspect=1
      endif
      weight=m*real(sd_n(n),RP)/area
      maxd=max(maxd,d)
      if(sd_rk(n)<VALID2INVALID)then
        sl_out=sl_out+weight
        cell=real(floor(sd_x(n)/sl_dx)*(JE-JS+1)+floor(sd_y(n)/sl_dx),RP)
        write(sl_au)sl_time,cell,real(sd_n(n),RP),m,d,rho,aspect,real(sdi%nmono(n),RP),real(sd_liqice(n),RP)
      else
        suspended=suspended+weight;live=live+1
      endif
    enddo
    sl_peak=max(sl_peak,real(live,RP));sl_maxd=max(sl_maxd,maxd)
    vapor=sum(dens(KS:KS+sl_layers-1,IS:IE,JS:JE)*qtrc(KS:KS+sl_layers-1,IS:IE,JS:JE,I_QV)) &
      *sl_depth/sl_layers/real((IE-IS+1)*(JE-JS+1),RP)
    source=sl_supply
    if(sl_mode==1)source=sl_dep
    res=suspended+sl_out+vapor-sl_vinitial-sl_in-source
    sl_maxerr=max(sl_maxerr,abs(res))
    if(mod(sl_calls,60)==0 .or. sl_time>=sl_end-1e-7_RP .or. sl_calls==1)then
      call sdm_rhot_qtrc2p_t(rhot,qtrc,dens,pres,temp)
      write(sl_lu,'(16(ES25.16,1X))')sl_time,sl_in,sl_dep,sl_supply,sl_out,suspended,vapor,res,sl_maxerr, &
        minval(temp(KS:KS+sl_layers-1,IS:IE,JS:JE)),maxval(temp(KS:KS+sl_layers-1,IS:IE,JS:JE)), &
        sl_maxd,real(live,RP),sl_peak,sl_removed_heat
      flush(sl_lu);flush(sl_au)
    endif
  end subroutine
