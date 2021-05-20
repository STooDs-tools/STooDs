module STOODS_mcmc
!~**********************************************************************
!~* Purpose: MCMC samplers for models in Space, Time Or Other
!~*          DimensionS (STOODS)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon, Univerity of Adelaide
!~**********************************************************************
!~* Created: 12/07/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~*     1. XXX
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. XXX
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
Public::OATsampler,Config_MCMC_Read

type::indices ! vector or indices
    integer(mik),allocatable::indx(:)
end type indices

type::vector ! vector or reals
    real(mrk),allocatable::v(:)
end type vector

type::matrix ! matrix or reals
    real(mrk),allocatable::m(:,:)
end type matrix

type::valScope ! scope of a process value
    type(indices),allocatable::val(:) ! size size(process%val); list of data affected by process%val(k)
end type valScope

! Object used to store relevant information during sampling
! (e.g.old/new values of inference functions, inverse covariances, etc.)
type::samplingInfo
    real(mrk)::fullPost                          ! Full posterior, only used with slow sampler
    real(mrk),allocatable::lkh(:)                ! size model%dat%n; lkh contribution of each data
    real(mrk),allocatable::prior(:)              ! size model%npar; prior for each parameter
    real(mrk),allocatable::H(:)                  ! size model%npro; H contribution of each process
    type(vector),allocatable::hyperPrior(:)      ! size model%npro; priors for hyperparameters of each process
    type(vector),allocatable::distancePrior(:)   ! size model%nDim; priors for distance parameters of each dimension
    type(valScope),allocatable::proScope(:)      ! size model%npro: scope of values of each process
end type samplingInfo
type(samplingInfo)::SAMI

! MCMC tunings
type,public:: MCMCTunings
    integer(mik)::Nsim=10000,Nadapt=10,StopAdapt=5000,Nburn=0,Nslim=1
    real(mrk):: minMR=0.1_mrk,maxMR=0.5_mrk,DownMult=0.9_mrk,UpMult=1.1_mrk
    logical::doSpeedUp=.true.,doJumpReport=.false.,doMoveRateReport=.false.
    character(250)::file='MCMC.txt',sep=';',jumpFile='',moveRateFile=''
end type MCMCTunings

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine OATsampler(model,mcmc,err,mess)
!^**********************************************************************
!^* Purpose: One-At-a-Time sampler
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.mcmc
!^* INOUT
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use STOODS_types, only:stoodsModel,unfold,unfoldNames
type(MCMCTunings),intent(in)::mcmc
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='OATsampler'
integer(mik)::i,j,unt,ncol,iter,unt_j,unt_mr
real(mrk)::vector(model%nVec),prior,lkh,H,post,denom
character(250)::head(model%nVec)

err=0;mess='';
ncol=model%nVec+4 ! add 4 columns for prior, lkh, H and post
!==============================================================
! INITIALISATION
call initSAMI(model,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
!==============================================================
! SAMPLING
write(*,*) '===================================='
write(*,*) 'Started sampling...'
write(*,*) '===================================='
! Initialise MCMC files
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(mcmc%file), status='replace', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
if(mcmc%doJumpReport) then
    call getSpareUnit(unt_j,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt_j,file=trim(mcmc%jumpFile), status='replace', iostat=err)
    if(err/=0) then;mess=trim(procname)//':open error';return;endif
endif
if(mcmc%doMoveRateReport) then
    call getSpareUnit(unt_mr,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt_mr,file=trim(mcmc%moveRateFile), status='replace', iostat=err)
    if(err/=0) then;mess=trim(procname)//':open error';return;endif
endif
! unfold and write headers
call unfoldNames(model,head,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
call writeOneLine_str(unt,(/head,'prior','lkh','H','post'/),mcmc%sep)
if(mcmc%doJumpReport) call writeOneLine_str(unt_j,head,mcmc%sep)
if(mcmc%doMoveRateReport) call writeOneLine_str(unt_mr,head,mcmc%sep)
! unfold and write starting value if needed
if(mcmc%Nburn==0) then
    call unfold(model,vector,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    call getPostComponents(prior,lkh,H,post)
    call writeOneLine(unt,(/vector,prior,lkh,H,post /),mcmc%sep)
endif
if(mcmc%doJumpReport) then
    call unfold(model,vector,err,mess,'jump')
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    call writeOneLine(unt_j,vector,mcmc%sep)
endif
do iter=1,mcmc%Nsim
    !==============================================================
    ! Update parameters
    do i=1,model%nPar
        call updatePar(model=model,i=i,doSpeedUp=mcmc%doSpeedUp,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    enddo
    !==============================================================
    ! Update hyperparameters and process values
    do i=1,model%nPro
         ! Update hyperparameters
        do j=1,model%process(i)%hyper%nPar
            call updateHyperPar(model=model,ihyper=j,ipro=i,doSpeedUp=mcmc%doSpeedUp,&
                                err=err,mess=mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        enddo
       ! Update process values
        do j=1,size(model%process(i)%val)
            call updateProcess(model=model,ival=j,ipro=i,doSpeedUp=mcmc%doSpeedUp,&
                               err=err,mess=mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        enddo
    enddo
    !==============================================================
    ! 2DO: Update distance parameters
    do i=1,model%nDim
        if(model%dims(i)%dist%nPar>0) then
            mess=trim(procname)//':distance parameters not supported yet'
            err=1;return
        endif
    enddo
    !==============================================================
    ! Save iteration
    if(iter>mcmc%Nburn .and. mod(iter-(mcmc%Nburn+1),mcmc%Nslim)==0) then
        ! unfold
        call unfold(model,vector,err,mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        call getPostComponents(prior,lkh,H,post)
        call writeOneLine(unt,(/vector,prior,lkh,H,post/),mcmc%sep)
    endif
    !==============================================================
    ! update jump sdev
    if(mod(iter,mcmc%Nadapt)==0) then
        if(mcmc%doMoveRateReport) then
            call unfold(model,vector,err,mess,'move')
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            if(iter<=mcmc%StopAdapt) then
                denom=real(mcmc%nAdapt,mrk)
            else
                denom=real(iter-mcmc%StopAdapt,mrk)
            endif
            call writeOneLine(unt_mr,vector/denom,mcmc%sep)
        endif
        write(*,*) 'Done:',100*real(iter)/real(mcmc%Nsim),'%'
        if(iter<=mcmc%StopAdapt) call updateJumpSD(mcmc,model)
        if(mcmc%doJumpReport) then
            call unfold(model,vector,err,mess,'jump')
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            call writeOneLine(unt_j,vector,mcmc%sep)
        endif
    endif
enddo
write(*,*) '===================================='
write(*,*) 'Finished sampling!'
write(*,*) '===================================='
close(unt)
if(mcmc%doJumpReport) close(unt_j)
if(mcmc%doMoveRateReport) close(unt_mr)
end subroutine OATsampler
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Config_MCMC_Read(workspace,file,mcmc,err,mess)
!^**********************************************************************
!^* Purpose: Read config file for MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 11/03/2020
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.workspace
!^*    2.file
!^* OUT
!^*    1.mcmc, MCMCTunings object
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::workspace,file
type(MCMCTunings),intent(out)::mcmc
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_MCMC_Read'
integer(mik)::unt
character(250)::cfile

err=0;mess=''
!==============================================================
! Open
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(workspace)//trim(file), status='old', iostat=err)
if(err/=0) then;mess=trim(procname)//':open error';return;endif
!==============================================================
! N's
read(unt,*,iostat=err) mcmc%Nsim
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%Nadapt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%StopAdapt
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%Nburn
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%Nslim
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Adaption tunings
read(unt,*,iostat=err) mcmc%minMR
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%maxMR
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%DownMult
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) mcmc%UpMult
if(err/=0) then;mess=trim(procname)//':read error';return;endif
!==============================================================
! Others
read(unt,*,iostat=err) cfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
mcmc%file=trim(workspace)//trim(cfile)
read(unt,*,iostat=err) mcmc%doSpeedUp
if(err/=0) then;mess=trim(procname)//':read error';return;endif
read(unt,*,iostat=err) cfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(trim(cfile)/='') then
    mcmc%doJumpReport=.true.
    mcmc%jumpFile=trim(workspace)//trim(cfile)
endif
read(unt,*,iostat=err) cfile
if(err/=0) then;mess=trim(procname)//':read error';return;endif
if(trim(cfile)/='') then
    mcmc%doMoveRateReport=.true.
    mcmc%moveRateFile=trim(workspace)//trim(cfile)
endif

end subroutine Config_MCMC_Read
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*********!
! PRIVATE !
!*********!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updatePar(model,i,doSpeedUp,err,mess)
!^**********************************************************************
!^* Purpose: Metropolis update for the ith parameter
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.i, parameter index
!^*    2.doSpeedUp? if .false., might be extremely slow!
!^* INOUT
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel,FIX_str
use STOODS_inference, only:getLogPrior1,getLogLkh1,getLogPost
use Distribution_tools,only:generate
integer(mik), intent(in)::i
logical,intent(in)::doSpeedUp
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='updatePar'
real(mrk)::currentVal,jump,lkh(model%dat%n),prior,logratio,u,oldPost,newPost
logical::feas,isnull
integer(mik)::j,n

err=0;mess='';
if(trim(model%par(i)%priorDist)==FIX_str) return ! fixed parameter

!==============================================================
! STEP 1: generate candidate
! remember current value
currentVal=model%par(i)%val
! generate jump
call Generate(DistId='Gaussian',par=(/0._mrk,model%par(i)%jump/),gen=jump,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) then;mess=trim(procname)//':jump generation unfeasible';return;endif
! update par value
model%par(i)%val=currentVal+jump
!==============================================================
! STEP 2: update prior and likelihood
if(doSpeedUp) then
    oldPost=SAMI%prior(i)+sum(SAMI%lkh)
    call getLogPrior1(par=model%par(i),logp=prior,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
        model%par(i)%val=currentVal
        return
    endif
    do j=1,model%dat%n
        call getLogLkh1(model=model,i=j,logp=lkh(j),feas=feas,isnull=isnull,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
            model%par(i)%val=currentVal
            return
        endif
    enddo
    newPost=prior+sum(lkh)
else
    oldPost=SAMI%fullPost
    call getLogPost(model=model,lpost=newPost,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
        model%par(i)%val=currentVal
        return
    endif
endif
!==============================================================
! STEP 3: apply Metropolis decision rule
logratio=newPost - oldPost
call Generate(DistId='Uniform',par=(/0._mrk,1._mrk/),gen=u,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(log(u)<=logratio) then ! accept and update everything
    model%par(i)%move=model%par(i)%move+1
    if(doSpeedUp) then
        SAMI%fullPost=SAMI%fullPost-SAMI%prior(i)-sum(SAMI%lkh)+prior+sum(lkh)
        SAMI%prior(i)=prior
        SAMI%lkh=lkh
    else
        SAMI%fullPost=newPost
    endif
else ! reject and revert
    model%par(i)%val=currentVal
endif
end subroutine updatePar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updateProcess(model,ival,ipro,doSpeedUp,err,mess)
!^**********************************************************************
!^* Purpose: Metropolis update for the (ival)th value of the (ipro)th process
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.ival, value index
!^*    2.ipro, process index
!^*    3.doSpeedUp? if .false., might be extremely slow!
!^* INOUT
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel,applyConstraints
use STOODS_inference, only:getLogLkh1,getLogH1,getLogPost
use Distribution_tools,only:generate
use MultivariateDistribution_tools,only:getLogPdf_singleCompUpdate
! arguments
integer(mik), intent(in)::ival,ipro
logical,intent(in)::doSpeedUp
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='updateProcess'
real(mrk)::currentVal,jump,oldH,newH,oldLkh,newLkh,logp,logratio,u,oldPost,newPost
real(mrk),allocatable::lkh(:)
integer(mik),allocatable::indx(:)
logical::feas,isnull
integer(mik)::i,n
real(mrk)::m,s2inv
real(mrk)::mu(size(model%process(ipro)%val)),&
           sk(size(model%process(ipro)%val)),&
           currentProc(size(model%process(ipro)%val)),&
           currentLogps(size(model%process(ipro)%val))

err=0;mess='';
! exit if values are constrained
if(ival > size(model%process(ipro)%val) - model%process(ipro)%constraint) return
!==============================================================
! STEP 1: generate candidate
! remember current values
currentProc=model%process(ipro)%val
currentVal=currentProc(ival)
if(allocated(model%process(ipro)%hyper%mdpar%logps)) then
    currentLogps=model%process(ipro)%hyper%mdpar%logps
endif
! generate jump
call Generate(DistId='Gaussian',par=(/0._mrk,model%process(ipro)%jump(ival)/),gen=jump,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) then;mess=trim(procname)//':jump generation unfeasible';return;endif
! update par value
model%process(ipro)%val(ival)=currentVal+jump
call applyConstraints(model%process(ipro),feas)
if(.not. feas) then ! revert process values and get out
   model%process(ipro)%val=currentProc
   return
endif
if(doSpeedUp) then
    !==============================================================
    ! STEP 2: update likelihood
    n=size(SAMI%proScope(ipro)%val(ival)%indx)
    if(allocated(lkh)) deallocate(lkh);allocate(lkh(n))
    if(allocated(indx)) deallocate(indx);allocate(indx(n))
    if(n>0) indx=SAMI%proScope(ipro)%val(ival)%indx
    do i=1,n
        call getLogLkh1(model=model,i=indx(i),logp=lkh(i),feas=feas,isnull=isnull,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) then ! revert process values and get out
           model%process(ipro)%val=currentProc
           return
        endif
    enddo
    if(n==0) then
        oldLkh=0._mrk;newLkh=0._mrk
    else
        oldLkh=sum(SAMI%lkh(indx));newLkh=sum(lkh)
    endif
    !==============================================================
    ! STEP 3: update hyperdist
    oldH=SAMI%H(ipro)
    call getLogPdf_singleCompUpdate(DistId=model%process(ipro)%hyper%name,&
                                    x=currentProc,&
                                    par=model%process(ipro)%hyper%mdpar,&
                                    oldLogp=oldH,iComp=ival,inc=jump,&
                                    newLogp=newH,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert process values and get out
        model%process(ipro)%val=currentProc
        if(allocated(model%process(ipro)%hyper%mdpar%logps)) then
            model%process(ipro)%hyper%mdpar%logps=currentLogps
        endif
        return
    endif
    oldPost=oldLkh+oldH
    newPost=newLkh+newH
else
    !==============================================================
    ! STEP 2&3: update full posterior
    oldPost=SAMI%fullPost
    call getLogPost(model=model,lpost=newPost,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert process values and get out
        model%process(ipro)%val=currentProc
        if(allocated(model%process(ipro)%hyper%mdpar%logps)) then
            model%process(ipro)%hyper%mdpar%logps=currentLogps
        endif
        return
    endif
endif
!==============================================================
! STEP 4: apply Metropolis decision rule
logratio=newPost-oldPost
call Generate(DistId='Uniform',par=(/0._mrk,1._mrk/),gen=u,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(log(u)<=logratio) then ! accept and update everything
    model%process(ipro)%move(ival)=model%process(ipro)%move(ival)+1
    if(doSpeedUp) then
        SAMI%fullPost=SAMI%fullPost-sum(SAMI%lkh(indx))-SAMI%H(ipro)+sum(lkh)+newH
        SAMI%lkh(indx)=lkh
        SAMI%H(ipro)=newH
    else
        SAMI%fullPost=newPost
    endif
else ! reject and revert
    model%process(ipro)%val=currentProc
    if(allocated(model%process(ipro)%hyper%mdpar%logps)) then
        model%process(ipro)%hyper%mdpar%logps=currentLogps
    endif
endif
end subroutine updateProcess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updateHyperPar(model,ihyper,ipro,doSpeedUp,err,mess)
!^**********************************************************************
!^* Purpose: Metropolis update for the (ihyper)th hyperpar of the (ipro)th process
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.ihyper, hyperparameter index
!^*    2.ipro, process index
!^*    3.doSpeedUp? if .false., might be extremely slow!
!^* INOUT
!^*    1.model, stoods model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel,FIX_str
use STOODS_inference, only:getLogPrior1,getLogH1,getLogPost
use Distribution_tools,only:generate
use MultivariateDistribution_tools,only:mDistPar
! arguments
integer(mik), intent(in)::ihyper,ipro
logical,intent(in)::doSpeedUp
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='updateHyperPar'
real(mrk)::currentVal,jump,oldPrior,newPrior,oldH,newH,logratio,u,oldPost,newPost
logical::feas,isnull
type(mDistPar)::currentPar

err=0;mess='';
if(trim(model%process(ipro)%hyper%par(ihyper)%priorDist)==FIX_str) return ! fixed hyperparameter

!==============================================================
! STEP 1: generate candidate
! remember current values
currentVal=model%process(ipro)%hyper%par(ihyper)%val
currentPar=model%process(ipro)%hyper%mdpar
! generate jump
call Generate(DistId='Gaussian',par=(/0._mrk,model%process(ipro)%hyper%par(ihyper)%jump/),&
              gen=jump,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) then;mess=trim(procname)//':jump generation unfeasible';return;endif
! update par value
model%process(ipro)%hyper%par(ihyper)%val=currentVal+jump
if(doSpeedUp) then
    oldPost=SAMI%hyperPrior(ipro)%v(ihyper)+SAMI%H(ipro)
    !==============================================================
    ! STEP 2: update prior and hierarchical part
    oldPrior=SAMI%hyperPrior(ipro)%v(ihyper)
    call getLogPrior1(par=model%process(ipro)%hyper%par(ihyper),logp=newPrior,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
        model%process(ipro)%hyper%par(ihyper)%val=currentVal
        return
    endif
    oldH=SAMI%H(ipro)
    call getLogH1(model=model,i=ipro,logp=newH,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
        model%process(ipro)%hyper%par(ihyper)%val=currentVal
        model%process(ipro)%hyper%mdpar=currentPar
        return
    endif
    newPost=newPrior+newH
else
    !==============================================================
    ! STEP 2: update full posterior
    oldPost=SAMI%fullPost
    call getLogPost(model=model,lpost=newPost,feas=feas,isnull=isnull,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then ! revert parameter value and get out
        model%process(ipro)%hyper%par(ihyper)%val=currentVal
        model%process(ipro)%hyper%mdpar=currentPar
        return
    endif
endif
!==============================================================
! STEP 3: apply Metropolis decision rule
logratio=newPost-oldPost
call Generate(DistId='Uniform',par=(/0._mrk,1._mrk/),gen=u,feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(log(u)<=logratio) then ! accept and update everything
    model%process(ipro)%hyper%par(ihyper)%move=model%process(ipro)%hyper%par(ihyper)%move+1
    if(doSpeedUp) then
        SAMI%fullPost=SAMI%fullPost-SAMI%hyperPrior(ipro)%v(ihyper)-SAMI%H(ipro)+newprior+newH
        SAMI%hyperPrior(ipro)%v(ihyper)=newprior
        SAMI%H(ipro)=newH
    else
        SAMI%fullPost=newPost
    endif
else ! reject and revert
    model%process(ipro)%hyper%par(ihyper)%val=currentVal
    model%process(ipro)%hyper%mdpar=currentPar
endif
end subroutine updateHyperPar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updateJumpSD(mcmc,model)
!^**********************************************************************
!^* Purpose: Update jump SD based on move rates
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 15/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.mcmc, mcmc tunings
!^* INOUT
!^*    1.model, stoods model object
!^**********************************************************************
use STOODS_types, only:stoodsModel,FIX_str
type(MCMCTunings),intent(in)::mcmc
type(stoodsModel),intent(inout)::model
! locals
integer(mik)::i,j,n

! Parameters
do i=1,model%nPar
    if(trim(model%par(i)%priorDist)==FIX_str) cycle ! fixed parameter
    call updateJump1(mcmc,model%par(i)%jump,model%par(i)%move)
enddo
! Processes
do i=1,model%nPro
    ! process values
    n=size(model%process(i)%val)
    do j=1,n
        if(j>n-model%process(i)%constraint) cycle ! constrained value
        call updateJump1(mcmc,model%process(i)%jump(j),model%process(i)%move(j))
    enddo
    ! hyperparameters
    n=model%process(i)%hyper%nPar
    do j=1,n
        if(trim(model%process(i)%hyper%par(j)%priorDist)==FIX_str) cycle ! fixed parameter
        call updateJump1(mcmc,model%process(i)%hyper%par(j)%jump,model%process(i)%hyper%par(j)%move)
    enddo
enddo
! Dimensions
do i=1,model%nDim
    n=model%dims(i)%dist%nPar
    do j=1,n
        if(trim(model%dims(i)%dist%par(j)%priorDist)==FIX_str) cycle ! fixed parameter
        call updateJump1(mcmc,model%dims(i)%dist%par(j)%jump,model%dims(i)%dist%par(j)%move)
    enddo
enddo

end subroutine updateJumpSD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine updateJump1(mcmc,jump,move)
!^**********************************************************************
!^* Purpose: Update a single jump SD based on number of moves
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 15/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.mcmc, mcmc tunings
!^* INOUT
!^*    1.jump, jump standard deviation
!^*    2.move, number of moves
!^**********************************************************************
type(MCMCTunings),intent(in)::mcmc
real(mrk),intent(inout)::jump,move
! locals
real(mrk)::MR

MR=move/real(mcmc%nAdapt,mrk)
if(MR<=mcmc%minMR) jump=jump*mcmc%DownMult
if(MR>=mcmc%maxMR) jump=jump*mcmc%UpMult
move=0._mrk

end subroutine updateJump1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine initSAMI(model,err,mess)
!^**********************************************************************
!^* Purpose: initialise SAMI object (sampling information)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1.model, stoods model object, to be updated
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:number_string
use MultivariateDistribution_tools,only:getLogPdf,mDistPar
use STOODS_types, only:stoodsModel
use STOODS_inference, only:getLogLkh1,getLogPrior1,getLogH1
use STOODS_formulas,only:getHyperPar
! arguments
type(stoodsModel),intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='initSAMI'
integer(mik)::i,j,k,n,m
logical::feas,isnull,singular
logical,allocatable::mask(:)
real(mrk)::fooSpar(0),signDet
!type(mDistPar)::gpar

err=0;mess='';

!==============================================================
! Allocate
if(allocated(SAMI%lkh)) deallocate(SAMI%lkh);allocate(SAMI%lkh(model%dat%n))
if(allocated(SAMI%prior)) deallocate(SAMI%prior);allocate(SAMI%prior(model%nPar))
if(allocated(SAMI%H)) deallocate(SAMI%H);allocate(SAMI%H(model%nPro))
if(allocated(SAMI%hyperPrior)) deallocate(SAMI%hyperPrior);allocate(SAMI%hyperPrior(model%nPro))
if(allocated(SAMI%distancePrior)) deallocate(SAMI%distancePrior);allocate(SAMI%distancePrior(model%nDim))
!==============================================================
! Get Likelihood contributions
do i=1,model%dat%n
    call getLogLkh1(model=model,i=i,logp=SAMI%lkh(i),feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then
        mess=trim(procname)//':'//'likelihood-unfeasible data:'//trim(number_string(i))
        err=1;return
    endif
enddo
!==============================================================
! Get Hierarchical contributions
do i=1,model%nPro
    call getLogH1(model=model,i=i,logp=SAMI%H(i),feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then
        mess=trim(procname)//':'//'Hierarchical-unfeasible process:'//trim(model%process(i)%name)
        err=1;return
    endif
enddo
!==============================================================
! Get Prior contributions
! parameter priors
do i=1,model%nPar
    call getLogPrior1(par=model%par(i),logp=SAMI%prior(i),feas=feas,isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) then
        mess=trim(procname)//':'//'Prior-unfeasible:'//trim(model%par(i)%name)
        err=1;return
    endif
enddo
! hyper-prior
do i=1,model%nPro
    m=model%process(i)%hyper%nPar
    if(allocated(SAMI%hyperPrior(i)%v)) deallocate(SAMI%hyperPrior(i)%v)
    allocate(SAMI%hyperPrior(i)%v(m))
    do j=1,m
        call getLogPrior1(par=model%process(i)%hyper%par(j),logp=SAMI%hyperPrior(i)%v(j),&
                          feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) then
            mess=trim(procname)//':'//'Prior-unfeasible:'//trim(model%process(i)%hyper%par(j)%name)
            err=1;return
        endif
    enddo
enddo
! distance-prior
do i=1,model%nDim
    m=model%dims(i)%dist%nPar
    if(allocated(SAMI%distancePrior(i)%v)) deallocate(SAMI%distancePrior(i)%v)
    allocate(SAMI%distancePrior(i)%v(m))
    do j=1,m
        call getLogPrior1(par=model%dims(i)%dist%par(j),logp=SAMI%distancePrior(i)%v(j),&
                          feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) then
            mess=trim(procname)//':'//'Prior-unfeasible:'//trim(model%dims(i)%dist%par(j)%name)
            err=1;return
        endif
    enddo
enddo
!==============================================================
! Full posterior
SAMI%fullPost=sum(SAMI%lkh) + sum(SAMI%H) + sum(SAMI%prior)
do i=1,size(SAMI%hyperPrior);SAMI%fullPost=SAMI%fullPost+sum(SAMI%hyperPrior(i)%v);enddo
do i=1,size(SAMI%distancePrior);SAMI%fullPost=SAMI%fullPost+sum(SAMI%distancePrior(i)%v);enddo
!==============================================================
! Work out scope of processes values
if(allocated(SAMI%proScope)) deallocate(SAMI%proScope);allocate(SAMI%proScope(model%nPro))
do i=1,model%nPro
    n=size(model%process(i)%val)
    k=model%proDimIndx(i) ! This process is associated with the kth dimension
    if(allocated(SAMI%proScope(i)%val)) deallocate(SAMI%proScope(i)%val)
    allocate(SAMI%proScope(i)%val(n))
    if(allocated(mask)) deallocate(mask);allocate(mask(model%dat%n))
    do j=1,n
        mask=(model%dat%idims(:,k)==j) ! data associated with the jth value of the process
        ! additional scope induced by contraints
        do m=1,model%process(i)%constraint
            mask= mask .or. ( model%dat%idims(:,k) == (n-m+1) )
        enddo
        if(allocated(SAMI%proScope(i)%val(j)%indx)) deallocate(SAMI%proScope(i)%val(j)%indx)
        allocate(SAMI%proScope(i)%val(j)%indx(count(mask)))
        SAMI%proScope(i)%val(j)%indx=pack((/(m,m=1,model%dat%n)/),mask)
    enddo
enddo
!==============================================================

end subroutine initSAMI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getPostComponents(prior,lkh,H,post)
!^**********************************************************************
!^* Purpose: get posterior components from SAMI object
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 15/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* OUT
!^*    1.prior: log-prior
!^*    2.lkh: log-lkh
!^*    3.H: log-hierarchical
!^*    4.post: log-post
!^**********************************************************************
real(mrk),intent(out)::prior,lkh,H,post
! locals
integer(mik)::i

lkh=sum(SAMI%lkh)
H=sum(SAMI%H)
prior=sum(SAMI%prior)
do i=1,size(SAMI%hyperPrior);prior=prior+sum(SAMI%hyperPrior(i)%v);enddo
do i=1,size(SAMI%distancePrior);prior=prior+sum(SAMI%distancePrior(i)%v);enddo
post=lkh+H+prior

end subroutine getPostComponents
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine writeOneLine(unt,vect,sep)
integer(mik),intent(in)::unt
real(mrk),intent(in)::vect(:)
character(*),intent(in)::sep
! locals
integer(mik)::i,n
character(250)::sep0

n=size(vect)
if(trim(sep)=="\t")then;sep0=achar(9);else;sep0=trim(sep);endif
do i=1,n
    if(i<n) then
        write(unt,'(e14.6,A)',advance='NO') vect(i),trim(sep0)
    else ! no separator for the last item
        write(unt,'(e14.6,A)',advance='NO') vect(i)
    endif
enddo
write(unt,'(A)',advance='YES') ""
end subroutine writeOneLine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine writeOneLine_str(unt,vect,sep)
integer(mik),intent(in)::unt
character(*),intent(in)::vect(:),sep
! locals
integer(mik)::i,n
character(250)::sep0

n=size(vect)
if(trim(sep)=="\t")then;sep0=achar(9);else;sep0=trim(sep);endif
do i=1,n
    if(i<n) then
        write(unt,'(A,A)',advance='NO') trim(vect(i)),trim(sep0)
    else ! no separator for the last item
        write(unt,'(A,A)',advance='NO') trim(vect(i))
    endif
enddo
write(unt,'(A)',advance='YES') ""
end subroutine writeOneLine_str
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module STOODS_mcmc
