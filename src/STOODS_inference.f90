module STOODS_inference
!~**********************************************************************
!~* Purpose: Inference functions for models in Space, Time Or Other
!~*          DimensionS (STOODS)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon, Univerity of Adelaide
!~**********************************************************************
!~* Created: 08/07/2019
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
Public::getLogLkh,getLogPrior,getLogH,getLogPost,&   ! full
        getLogLkh1,getLogPrior1,getLogH1             ! one data / parameter / process

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogLkh(model,ll,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute full log-likelihood
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 08/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Optimize computing time by avoiding potentially useless
!^*           applications of formulas.
!^*           Impovement might be negligible given other computational
!^*           bottlenecks though.
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^* OUT
!^*    1.ll, log-likelihood
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel
type(stoodsModel),intent(in)::model
real(mrk), intent(out)::ll
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogLkh'
integer(mik)::i
real(mrk)::logp

err=0;mess='';feas=.true.;isnull=.false.;ll=0._mrk

do i=1,model%dat%n ! 2DO: smarter looping to avoid useless applications of formulas
    call getLogLkh1(model,i,logp,feas,isnull,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) return
    ll=ll+logp
enddo

end subroutine getLogLkh
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogPrior(model,lp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute full log-prior
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^* OUT
!^*    1.lp, prior log-pdf
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel
type(stoodsModel),intent(in)::model
real(mrk), intent(out)::lp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogPrior'
integer(mik)::i,j
real(mrk)::logp

err=0;mess='';feas=.true.;isnull=.false.;lp=0._mrk

! parameters
do i=1,model%nPar
    call getLogPrior1(model%par(i),logp,feas,isnull,err,mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) return
    lp=lp+logp
enddo

! hyper-parameters
do i=1,model%nPro
    do j=1,model%process(i)%hyper%nPar
        call getLogPrior1(model%process(i)%hyper%par(j),logp,feas,isnull,err,mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) return
        lp=lp+logp
    enddo
enddo

! dimensions: distance parameters
do i=1,model%nDim
    do j=1,model%dims(i)%dist%nPar
        call getLogPrior1(model%dims(i)%dist%par(j),logp,feas,isnull,err,mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull ) return
        lp=lp+logp
    enddo
enddo

end subroutine getLogPrior
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogH(model,lh,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute log of the hierarchical component
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 10/07/2019
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
!^*    1.lh, hierarchical log-pdf
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel
type(stoodsModel),intent(inout)::model
real(mrk), intent(out)::lh
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogH'
integer(mik)::i
real(mrk)::logp

err=0;mess='';feas=.true.;isnull=.false.;lh=0._mrk

do i=1,model%nPro
    call getLogH1(model,i,logp,feas,isnull,err,mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( (.not. feas) .or. isnull ) return
    lh=lh+logp
enddo

end subroutine getLogH
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogPost(model,lpost,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute full log-posterior
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 10/07/2019
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
!^*    1.lpost, posterior log-pdf
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel
type(stoodsModel),intent(inout)::model
real(mrk), intent(out)::lpost
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogPost'
real(mrk)::ll,lp,lh

err=0;mess='';feas=.true.;isnull=.false.;lpost=0._mrk
! Prior
call getLogPrior(model,lp,feas,isnull,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if( (.not. feas) .or. isnull ) return
lpost=lpost+lp
! Likelihood
call getLogLkh(model,ll,feas,isnull,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if( (.not. feas) .or. isnull ) return
lpost=lpost+ll
! hierarchical
call getLogH(model,lh,feas,isnull,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if( (.not. feas) .or. isnull ) return
lpost=lpost+lh

end subroutine getLogPost
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogLkh1(model,i,logp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute log-likelihood contribution for ith data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 12/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Optimize computing time by avoiding potentially useless
!^*           applications of formulas.
!^*           Impovement might be negligible given other computational
!^*           bottlenecks though.
!^**********************************************************************
!^* IN
!^*    1.model, stoods model object
!^* OUT
!^*    1.ll, log-likelihood
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel
use STOODS_formulas,only:getParentPar
use Distribution_tools,only:GetPdf,GetCdf
type(stoodsModel),intent(in)::model
integer(mik), intent(in)::i
real(mrk), intent(out)::logp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogLkh1'
integer(mik)::ivar,k
real(mrk)::par(sum(model%nParentPar)) ! not that depending on the variable, only the first model%nParentPar(ivar) should be used
character(250)::dist
real(mrk)::lp1,lp2

err=0;mess='';feas=.true.;isnull=.false.;logp=undefRN

par=undefRN
! retrieve parameters
call getParentPar(model=model,i=i,parentPar=par,&
                  feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) return
! compute log-pdf of obs
ivar=model%dat%ivar(i) ! Which variable?
k=model%nParentPar(ivar) ! Only use par(1:k)
dist=model%parentDist(ivar) ! Which distribution?
if(model%dat%cType(i)==0) then ! no or interval censoring
    if(model%dat%cWidth(i)==0._mrk) then ! no censoring
        call GetPdf(DistId=trim(dist),x=model%dat%y(i),par=par(1:k),&
                    loga=.true.,pdf=logp,feas=feas,isnull=isnull,&
                    err=err,mess=mess)
        return
    else ! interval censoring
        call GetCdf(DistId=trim(dist),x=model%dat%y(i)+model%dat%cWidth(i),&
                    par=par(1:k),cdf=lp2,feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) return
        call GetCdf(DistId=trim(dist),x=model%dat%y(i)-model%dat%cWidth(i),&
                    par=par(1:k),cdf=lp1,feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) return
        if(abs(lp2-lp1)==0._mrk) then
            isnull=.true.;return
        else
            logp=log(abs(lp2-lp1));return
        endif
    endif
endif
if(model%dat%cType(i)<0)then ! "less than" censoring
        call GetCdf(DistId=trim(dist),x=model%dat%y(i),&
                    par=par(1:k),cdf=lp2,feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) return
        if(lp2<=0._mrk) then
            isnull=.true.;return
        else
            logp=log(lp2);return
        endif
endif
if(model%dat%cType(i)>0)then ! "more than" censoring
        call GetCdf(DistId=trim(dist),x=model%dat%y(i),&
                    par=par(1:k),cdf=lp1,feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) return
        if(lp1>=1._mrk) then
            isnull=.true.;return
        else
            logp=log(1._mrk-lp1);return
        endif
endif
end subroutine getLogLkh1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogPrior1(par,logp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute log-prior for a single parameter
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.par, param object
!^* OUT
!^*    1.logp, prior log-pdf
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsParam,FIX_str
use Distribution_tools,only:GetPdf
type(stoodsParam),intent(in)::par
real(mrk), intent(out)::logp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogPrior1'

err=0;mess='';feas=.true.;isnull=.false.;logp=undefRN
if(trim(par%priorDist)==FIX_str) then
    logp=0._mrk
else
    call GetPdf(DistId=trim(par%priorDist),x=par%val,par=par%priorPar,&
            loga=.true.,pdf=logp,feas=feas,isnull=isnull,err=err,mess=mess)
endif

! Kill warnings before returning
if(err<0) err=0

end subroutine getLogPrior1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine getLogH1(model,i,logp,feas,isnull,err,mess)
!^**********************************************************************
!^* Purpose: Compute log of the hierarchical component for a single process
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon / University of Adelaide
!^**********************************************************************
!^* Created: 09/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1.model, stoods model object, to be updated
!^* IN
!^*    1.i, process index
!^* OUT
!^*    1.logp, log-prior pdf
!^*    2.feas, feasible?
!^*    3.isnull, is lkh=0?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use STOODS_types, only:stoodsModel,stoodsProcess
use STOODS_formulas,only:getHyperPar
use MultivariateDistribution_tools,only:GetLogPdf,mDistPar
! arguments
type(stoodsModel),intent(inout)::model
integer(mik),intent(in)::i
real(mrk), intent(out)::logp
logical, intent(out)::feas,isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='getLogH1'
type(stoodsProcess)::pro

err=0;mess='';feas=.true.;isnull=.false.;logp=undefRN
pro=model%process(i)
! retrieve parameters
call getHyperPar(model=model,i=i,&
                 feas=feas,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) return

! Compute pdf
call GetLogPdf(DistId=pro%hyper%name,x=pro%val,par=model%process(i)%hyper%mdpar,&
               logp=logp,feas=feas,isnull=isnull,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) return

end subroutine getLogH1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module STOODS_inference
