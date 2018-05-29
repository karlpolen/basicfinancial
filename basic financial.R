#loan formulas
loanbal=function(bal,i,pmt,n) {
  (bal*((1+i)^n))-((pmt/i)*(((1+i)^n)-1))
}
lpmt=function(bal,i,n) {
  (i*bal)/(1-((1+i)^-n))
}
intportion=function(bal,i,pmt,n) {
  i*loanbal(bal,i,pmt,n-1)
}
prinportion=function(bal,i,pmt,n) {
  pmt-intportion(bal,i,pmt,n)
}
loanamort=function(r=NULL,bal0=NULL,pmt=NULL,n=NULL,apr=FALSE,start=NULL,freq=1) {
  ans=list()
  risnull=is.null(r)
  bal0isnull=is.null(bal0)
  if(!bal0isnull) {
    if(is.zoo(bal0)) start=time(bal0)
    bal0=as.numeric(bal0)
  }
  pmtisnull=is.null(pmt)
  nisnull=is.null(n)
  if(1<sum(c(risnull,bal0isnull,pmtisnull,nisnull))) stop('loanamort error -- need to provide at least three parameters')
  n.f=n
  if(apr) n.f=n*freq
  if(!risnull) {
    if(apr) {
      r.f=r/freq
    } else {
      r.f=-1+(1+r)^(1/freq)
    }
  } else {
    cf=c(-bal0,rep(pmt,n.f))
    if(0<=sum(cf)) {
      rootrange=c(0,1.01) } else {
        rootrange=c(1,1000)
      }
    d=(uniroot(function(d) {sum(cf*d^(0:n.f))},rootrange))$root
    r.f=(1/d)-1
  }
  d=1/(1+r.f)
  f=1+r.f
  if(pmtisnull) pmt=(bal0*r.f)/(1-d^n.f)
  perp=pmt/r.f
  if(bal0isnull) bal0=perp-perp*(d^n)
  if(pmt<=(r.f*bal0)) stop(paste(pmt,r.f*bal0,'payment must be greater than interest'))
  if(nisnull) n.f= ceiling(log((1-(bal0*r.f)/pmt))/log(d))
  i=1:n.f
  bal=pmax(0,((bal0*f^i)-(((pmt*f^i)-pmt)/r.f)))
  balall=c(bal0,bal)
  int=balall[i]*r.f
  prin=-diff(balall)
  if(!is.null(start)) {
    bal=zooreg(bal,start=start+1/freq,freq=freq)
    int=zooreg(int,start=start+1/freq,freq=freq)
    prin=zooreg(prin,start=start+1/freq,freq=freq)
  }
  if(apr) {
    ans$r=r.f*freq
    ans$n=n.f/freq
  } else {
    ans$r=-1+((1+r.f)^freq)
    ans$n=n.f
  }
  ans$pmt=pmt
  ans$bal0=bal0
  ans$freq=freq
  ans$start=start
  ans$apr=apr
  ans$bal=bal
  ans$prin=prin
  ans$int=int
  return(ans)
}
#fv.z
fv.z=function(i,cf,freq=1,apr=FALSE,now=NULL,drop.after.now=TRUE,adjdisc=NULL) {
  if(!is.null(adjdisc)) apr=!adjdisc
  if(!is.zoo(cf)) {
    timeline=(1:length(cf))/freq
  } else {
    timeline=time(cf)
  }
  iszooreg=("Date"!=class(timeline))
  if(!iszooreg) {
    freq=365
  } else if ("yearmon" == class(timeline)) {
    freq=12
  } else if ("yearqtr" == class(timeline)) {
    freq=4
  }
  if (iszooreg&apr) {
    i=-1+(1+i/freq)^freq
  } else if (iszooreg&(!apr)) {
    i=i
  } else if ((!iszooreg)&apr) {
    i=i/freq
  } else if((!iszooreg)&(!apr)) {
    i=-1+((1+i)^(1/freq))
  }
  if (is.null(now)) now=timeline[length(timeline)] 
  if (drop.after.now) {
    cf=cf[timeline<=now]
    timeline=timeline[timeline<=now]
  }
  ind=as.numeric(now-timeline)
  d=(1+i)^ind
  sum(cf*d)
}

#npv
npv.z=function(i,cf,freq=1,apr=FALSE,now=NULL,drop.bef.now=TRUE,adjdisc=NULL) {
  if(!is.null(adjdisc)) apr=!adjdisc
  if(!is.zoo(cf)) {
    timeline=(1:length(cf))/freq
  } else {
    timeline=time(cf)
  }
  iszooreg=("Date"!=class(timeline))
  if(!iszooreg) {
    freq=365
  } else if ("yearmon" == class(timeline)) {
    freq=12
  } else if ("yearqtr" == class(timeline)) {
    freq=4
  }
  if (iszooreg&apr) {
    i=-1+(1+i/freq)^freq
  } else if (iszooreg&(!apr)) {
    i=i
  } else if ((!iszooreg)&apr) {
    i=i/freq
  } else if((!iszooreg)&(!apr)) {
    i=-1+((1+i)^(1/freq))
  }
  if (is.null(now)) now=timeline[1] 
  if (drop.bef.now) {
    cf=cf[timeline>=now]
    timeline=timeline[timeline>=now]
  }
  ind=as.numeric(timeline-now)
  d=(1+i)^ind
  sum(cf/d)
}
simple.interest.z=function(i,cf,start=time(cf)[1],end=lastinvec(time(cf))) {
  cf.h=cf[time(cf)>=start]
  cf.h=cf.h[time(cf.h)<=end]
  days=as.numeric(end-time(cf.h))
  int=i*(days/365)*cf.h
  return(sum(int))
}


# old.npv.z=function(i,cf.z,adjdisc=TRUE) {
#   freq=365
#   if(adjdisc) i=freq*(((1+i)^(1/freq))-1)
#   tdif=as.numeric(index(cf.z)-(index(cf.z)[1]))
#   d=(1+(i/freq))^tdif
#   sum(cf.z/d)
# }
# old.npv.znoadjust=function(i,cf.z) {
#   freq=365
#   if(!is.zoo(cf.z)) {warning("cash flow must be zoo object"); return(NA)}
#   if("Date"!=class(time(cf.z))) {warning("need Date class for zoo index"); return(NA)}
#   tdif=as.numeric(index(cf.z)-(index(cf.z)[1]))
#   d=(1+(i/freq))^tdif
#   sum(cf.z/d)
# }
irrgap.z=function(i,cf.z,end=NA) {
  freq=365
  if(!is.zoo(cf.z)) {warning("cash flow must be zoo object"); return(NA)}
  if("Date"!=class(time(cf.z))) {warning("need Date class for zoo index"); return(NA)}
  i=freq*(((1+i)^(1/freq))-1)
  if(!is.na(end)) cf.z=mergesum.z(cf.z,zoo(0,end))
  tdif=as.numeric(index(cf.z)-(index(cf.z)[1]))
  d=(1+(i/freq))^tdif
  v=sum(cf.z/d)
  -v*tail(d,1)
}
  
#irr
irr <- function(cf,rejectlm1=TRUE,frequ=1) { 
  res=polyroot(cf)
  res=Re(res[abs(Im(res))< 1e-6])
  res=1/res-1
  if(rejectlm1) res=res[res>=-1]
  if(0==length(res)){
    res=NA
    } else {
    res=sort(res)
    res=((1+res)^frequ)-1
  }
  return(res)
}
irr.z=function(cf.z,gips=FALSE) {
  irr.freq=1
  if(any(is.na(cf.z))) return(NA)
  if(sum(cf.z)==0) return (0)
  if(length(cf.z)<=1) return(NA)
  if(all(cf.z<=0)) return(NA)
  if(all(cf.z>=0)) return(NA)
  if(!is.zoo(cf.z)) {
    timediff=-1+1:length(cf.z)
  } else {
    timeline=time(cf.z)
    timediff=as.numeric(timeline-timeline[1])
    if ("Date"== class(timeline)) irr.freq=365
  }
  failirr=FALSE
  samesign=TRUE
  if (sum(cf.z)<0) {
    rangehi=.01
    rangelo=0
    i=0
     # low range on search for negative IRR is -100%
    while(i<1002&samesign&!failirr) {
      rangehi=rangelo
      rangelo=max(-10,rangelo-.01)
      hinpv=npv.znoadjust(rangehi,cf.z,irr.freq,timediff)
      lownpv=npv.znoadjust(rangelo,cf.z,irr.freq,timediff)
      samesign=sign(hinpv)==sign(lownpv)
      failirr=(is.nan(hinpv)|is.nan(lownpv)|is.na(samesign))
      i=i+1
    }} else {
      rangehi=0
      rangelo=-.01
      i=0
      # while hi range on search for positive IRR is 100,000%
      while(i<100000&samesign&!failirr) {
        rangelo=rangehi
        rangehi=rangehi+.01
        hinpv=npv.znoadjust(rangehi,cf.z,irr.freq,timediff)
        lownpv=npv.znoadjust(rangelo,cf.z,irr.freq,timediff)
        samesign=sign(hinpv)==sign(lownpv)
        failirr=(is.nan(hinpv)|is.nan(lownpv)|is.na(samesign))
        i=i+1
      }}
  if(failirr) return(NA)  
  npv1=lownpv
  npv2=hinpv
  if (sign(npv1)==sign(npv2)) return(NA)
  cf.n=as.numeric(cf.z)
  #calculate with uniroot if cash flow starts negative and ends positive otherwise do your own search
  if((cf.n[1]<0)&(cf.n[length(cf.n)]>0)) {
    ans=uniroot(npv.znoadjust,c(rangelo,rangehi),cf=cf.z,freq=irr.freq,tdiff=timediff)
    apr=ans$root } else {
      int1=rangelo
      int2=rangehi
      for (i in 1:40) {
        inta=mean(c(int1,int2))
        npva=npv.znoadjust(inta,cf.z,irr.freq,timediff)
        if(sign(npva)==sign(npv1)) {
          int1=inta
          npv1=npva
        } else {
          int2=inta
          npv2=npva
        }}
      apr=mean(int1,int2)  
    }
  # convert IRR to compounding at irr.freq interval
  ans=((1+(apr/irr.freq))^irr.freq)-1
  # convert IRR to GIPS compliant if requested
  if (gips) {
    if(cf.z[1]==0)  cf.z=cf.z[-1]
    dur=index(cf.z)[length(cf.z)]-index(cf.z)[1]
    if(dur<irr.freq) ans=(1+ans)^((as.numeric(dur))/irr.freq)-1
  }
  return (ans)
}
npv.znoadjust=function(i,cf.z,freq,tdiff) {
  d=(1+(i/freq))^tdiff
  sum(cf.z/d)
}

# irr.zold=function(cf.z,irr.freq=365,gips=FALSE) {
#   #if(any(cf.z==0)) cf.z=cf.z[cf.z!=0]
#   if(any(is.na(cf.z))) return(NA)
#   if(length(cf.z)<=1) return(NA)
#   if(all(cf.z<=0)) return(NA)
#   if(all(cf.z>=0)) return(NA)
#   if(sum(cf.z)==0) return (0)
#   #timespan=as.numeric(time(lastinvec(cf.z))-time(cf.z[1]))/irr.freq
#   neglower=-10
#   #if(timespan<1) neglower=-10
#   if (sum(cf.z)<=0) {
#     interval=c(neglower,0)} else {
#       interval=c(0,10)
#     }
#   npv1=npv.znoadjust(interval[1],cf.z)
#   npv2=npv.znoadjust(interval[2],cf.z)
#   if (sign(npv1)==sign(npv2)) return(NA)
#   cf.n=as.numeric(cf.z)
#   if((cf.n[1]<0)&(cf.n[length(cf.n)]>0)) {
#     ans=uniroot(npv.znoadjust,interval,cf=cf.z)
#     apr=ans$root } else {
#       int1=interval[1]
#       int2=interval[2]
#       npv1=npv.znoadjust(int1,cf.z)
#       npv2=npv.znoadjust(int2,cf.z)
#       for (i in 1:40) {
#         inta=mean(c(int1,int2))
#       npva=npv.znoadjust(inta,cf.z)
#       if(sign(npva)==sign(npv1)) {
#         int1=inta
#         npv1=npva
#       } else {
#         int2=inta
#         npv2=npva
#       }}
#     apr=mean(int1,int2)  
#     }
#   ans=((1+(apr/irr.freq))^irr.freq)-1
#   if (gips) {
#     if(cf.z[1]==0)  cf.z=cf.z[-1]
#     dur=lastinvec(index(cf.z))-index(cf.z)[1]
#     if(dur<irr.freq) ans=(1+ans)^((as.numeric(dur))/irr.freq)-1
#   }
#   return (ans)
# }
# old2.irr.z=function(cf.z,gips=FALSE) {
#   irr.freq=365
#   if(!is.zoo(cf.z)) {warning("cash flow must be zoo object"); return(NA)}
#   if("Date"!=class(time(cf.z))) {warning("need Date class for zoo index"); return(NA)}
#   if(any(is.na(cf.z))) return(NA)
#   if(length(cf.z)<=1) return(NA)
#   if(all(cf.z<=0)) return(NA)
#   if(all(cf.z>=0)) return(NA)
#   if(sum(cf.z)==0) return (0)
#   if (sum(cf.z)<0) {
#     rangehi=0
#     rangelo=-.01
#     i=0
#     while(i<10000&(sign(npv.znoadjust(rangehi,cf.z))==sign(npv.znoadjust(rangelo,cf.z)))) {
#       rangehi=rangelo
#       rangelo=rangelo-.01
#       i=i+1
#     }} else {
#       rangehi=.01
#       rangelo=0
#       i=0
#       while(i<10000&(sign(npv.znoadjust(rangehi,cf.z))==sign(npv.znoadjust(rangelo,cf.z)))) {
#         rangelo=rangehi
#         rangehi=rangehi+.01
#         i=i+1
#       }}
#   npv1=npv.znoadjust(rangelo,cf.z)
#   npv2=npv.znoadjust(rangehi,cf.z)
#   if (sign(npv1)==sign(npv2)) return(NA)
#   cf.n=as.numeric(cf.z)
#   if((cf.n[1]<0)&(cf.n[length(cf.n)]>0)) {
#     ans=uniroot(npv.znoadjust,c(rangelo,rangehi),cf=cf.z)
#     apr=ans$root } else {
#       int1=rangelo
#       int2=rangehi
#       for (i in 1:40) {
#         inta=mean(c(int1,int2))
#         npva=npv.znoadjust(inta,cf.z)
#         if(sign(npva)==sign(npv1)) {
#           int1=inta
#           npv1=npva
#         } else {
#           int2=inta
#           npv2=npva
#         }}
#       apr=mean(int1,int2)  
#     }
#   ans=((1+(apr/irr.freq))^irr.freq)-1
#   if (gips) {
#     if(cf.z[1]==0)  cf.z=cf.z[-1]
#     dur=lastinvec(index(cf.z))-index(cf.z)[1]
#     if(dur<irr.freq) ans=(1+ans)^((as.numeric(dur))/irr.freq)-1
#   }
#   return (ans)
# }


irr.z365=function(cf.z) {
  beg=min(index(cf.z))
  end=max(index(cf.z))
  cf.i=zooreg(rep(0,as.integer(1+end-beg)),start=beg,end=end)
  cf.mat=merge(cf.z,cf.i)
  cf=as.numeric(natozero(cf.mat$cf.z))
  res=irr(cf,frequ=365)
  return(res)
}
#tvpi
tvpi <- function(cf) {
  if(any(is.na(cf))) return(NA)
  if(0==sum(cf[cf<0])) {return(NA)}
  -sum(cf[cf>0])/sum(cf[cf<0])
}
#payback
payback <- function(cf) {
  cumcf=cumsum(cf)
  cumcfneg=rev(cumcf<0)
  pb=sum(cummax(cumcfneg))
  ans=pb
  if (pb==length(cf)) {ans=NA}
  ans
}
#performance report
invest.measures=function(cf,disc,f=Frequency) {
  cf.irr=100*irr(cf,frequ=f)
  name1="IRR"
  cf.npv=npv(disc/f,cf)
  name2= paste("NPV at",as.character(disc*100),"%")
  cf.tvpi=tvpi(cf)
  name3="Total Value Multiple to Invested"
  cf.payback=payback(cf)
  pbnames=c("years","2","3","Quarters","5","6","7","8","9","10","11","Months")
  name4=paste("Payback in",pbnames[f])
  ans=matrix(c(cf.irr,cf.npv,cf.tvpi,cf.payback),4,1)
  rownames(ans)=c(name1,name2,name3,name4)
  colnames(ans)=" "
  round(ans,2)
} 


# simple analysis of incentive pay in a multi-tier waterfall structure
incent=function(deal,return){
  if (is.null(dim(deal))) deal=matrix(deal,1,3,byrow=TRUE)
  colnames(deal)<-c("pref","carry","catchup")
  base=excess=matrix(0,0,length(return))
  incentfee=matrix(0,2,length(return))
  pref3=pref2=c(deal[,"pref"],10000)
  prefdif=diff(pref2)
  pref3[1]=0
  for (j in 1:dim(deal)[1]) {
    excess=rbind(excess,(pmax(0,pmin(prefdif[j],return-deal[j,"pref"]))))
    base=rbind(base,(pmin(pref2[j+1],return))*(return>pref3[j]))
    prior=(base[j,]>0)*colSums(incentfee)
    x=pmin(excess[j,]*deal[j,"catchup"],base[j,]*deal[j,"carry"]-prior)
    incentfee=rbind(incentfee,x)
  }
  colSums(incentfee)
}
gpcomp = function(dmat, ret, capital = 100, invcost = 100) {
  am = dmat$am
  pref = dmat$pref
  catchup = dmat$catchup
  carry = dmat$carry
  if (any(1 < c(catchup, carry))) 
    stop("catchup and carry must be stated as decimals<1")
  pref = c(pref, 1e+16)
  am = c(am, 0)
  stack = vector()
  lpcut = vector()
  typ = vector()
  nlayer = nrow(dmat)
  if (am[1] > 0) {
    stack = c(stack, am[1])
    lpcut = c(lpcut, 0)
    typ = c(typ, paste("Asset mgmt", 0))
  }
  if (capital > 0) {
    stack = c(stack, capital)
    lpcut = c(lpcut, 1)
    typ = c(typ, paste("Return of Capital"))
  }
  if (pref[1] > 0) {
    stack = c(stack, pref[1])
    lpcut = c(lpcut, 1)
    typ = c(typ, paste("Preferred Return", 1))
  }
  for (j in 1:nlayer) {
    if (am[j + 1] > 0) {
      stack = c(stack, am[j + 1])
      lpcut = c(lpcut, 0)
      typ = c(typ, paste("Asset Mgmt", j))
    }
    nextpref = pref[j + 1]
    lpsofar = sum(stack * lpcut) - capital
    lpshort = nextpref - lpsofar
    cu = catchup[j]
    cy = carry[j]
    catchuplayer = 0
    if (cu > cy) {
      catchuplayer = (lpsofar * cy)/(cu - cy)
      if (cu < 1) 
        catchuplayer = min(catchuplayer, lpshort/(1 - cu))
      stack = c(stack, catchuplayer)
      lpcut = c(lpcut, (1 - cu))
      typ = c(typ, paste("Catchup", j))
    }
    lpsofar = sum(stack * lpcut) - capital
    lpshort = nextpref - lpsofar
    carrylayer = lpshort/(1 - cy)
    if (carrylayer > 0) {
      stack = c(stack, carrylayer)
      lpcut = c(lpcut, (1 - cy))
      typ = c(typ, paste("Carry", j))
    }
  }
  ansmat = matrix(0, nrow = length(stack), ncol = length(ret))
  for (i in 1:length(ret)) {
    ansmat[, i] = wf(stack, ret[i])[-(1 + length(stack))]
  }
  ans = list()
  ans$lpshare = matrix(lpcut, nrow = length(stack), ncol = length(ret)) * 
    ansmat
  rownames(ans$lpshare) = typ
  ans$gpshare = ansmat - ans$lpshare
  rownames(ans$gpshare) = typ
  ans$grossreturn = 100 * (ret - invcost)/invcost
  ans$netreturn = 100 * (colSums(ans$lpshare) - capital)/capital
  ans$stack = stack
  ans$lpcut = lpcut
  return(ans)
}
# given waterfall w in dollars and available cash c, distribute the cash to
# the waterfall
wf = function(w, c) {
  x = c - cumsum(w)
  x[x > 0] = 0
  x = x + w
  x[x < 0] = 0
  c(x, c - sum(x))
}
testans = function(ans) {
  (colSums(ans$lpshare)) + (colSums(ans$gpshare))
}
#given waterfall w in dollars and available cash c, distribute the cash to the waterfall
waterfall=function(w,c){
  x=c-cumsum(w)
  x[x>0]=0
  x=x+w
  x[x<0]=0
  c(x,c-sum(x))
}

#constant growth
grow = function (c=1,t=term,g) {
  c*((1+g)^(0:(t-1)))
}
#last in vector
lastinvec=function(x,naremove=FALSE) {
  ans=tail(x,1)
  if(length(ans)> 0){
    if(naremove & is.na(ans)) ans=0
  }
  return(ans)
}
#first in vector
firstinvec=function(x) {x[1]}
#convert quarterly or monthly zoo object to annual zoo object
toannual=function(z,fun=sum) {
  ans=aggregate(z,as.character(floor(index(z))),fun)
  rownames(ans)=unique(as.character(floor(index(z))))
  ans
}
toannual.d=function(z,fun=sum,fiscal=.5) {
  ans=aggregate(z,as.yearqtr,fun)
  ans=aggregate(ans,as.character(floor(fiscal+index(ans))),fun)
  return(ans)
}
#convert zoo data to quarterly
toquarterly=function(z,fun=sum) {
  ans=aggregate(z,as.yearqtr(index(z)),fun)
  rownames(ans)=unique(as.character(as.yearqtr(index(z))))
  ans
}
#na to zero
natozero=function(x) {
  x[is.na(x)]=0
  x
}
#maketable is merge and change NA's to zero, establish names for variables
maketable=function(x,y,...,names=NULL) {
  ans=merge(x,y,...)
  ind=index(ans)
  ans=apply(ans,2,natozero)
  ans=zoo(ans,ind)
  if(!is.null(names)) colnames(ans)=names
  ans   
}
#monthdiff calculates difference in months between two dates represented c(year,mon)
monthdiff=function(b,e) {
  ans=12*(e[1]-b[1])
  ans+1-b[2]+e[2]
}
#monthadd calculates a new c(year,mon) from start and n=months to add
monthadd=function(start,n) {
  ans=start
  ans[2]=ans[2]+n-1
  ans[1]=ans[1]+floor(ans[2]/12)
  ans[2]=1+(ans[2]%%12)
  ans
}

mergesum.z=function(x,...) {
  mt=merge(x,...,fill=0)
  if(is.null(dim(mt))) return(mt)
  zoo::zoo(rowSums(mt),time(mt))
}

mergelist.z=function(z.list) {
  require(zoo)
  listlen=length(z.list)
  if (listlen==0) return(z.list)
  if (listlen==1) {
    ans=matrix(unlist(z.list),ncol=1)
    colnames(ans)=names(z.list)
    return(ans)
    }
  len.each=lapply(z.list,length)
  len.idx=which(len.each==0)
  for (i in len.idx) {
    z.list[[i]]=zoo(0,today())
  }
  ans=merge(z.list[[1]],z.list[[2]],fill=0)
  if (listlen==2) {
    colnames(ans)=names(z.list)
    return(ans)
  }
  for (k in 3:(listlen)) {
    ans=merge(ans,z.list[[k]],fill=0)
  }
  colnames(ans)=names(z.list)
  return(ans)
}
zoosum = function(...) {
  sumz = merge(...)
  if(is.null(dim(sumz))) return(sumz)
  zoo(rowSums(sumz, na.rm = TRUE), time(sumz))
}
addtotal.list = function(x) {
  total = do.call(zoosum, x)
  x$Total = total
  return(x)
}
merge0 = function(...) {
  merge(..., fill = 0)
}
#create a zoo object with global variable defaults
fullzoo=function(x=rep(0,Term),start=Start,end=End,frequency=Frequency,growth=0,lag=0,warn=FALSE) {
  if (length(x)==1) {
    vals=c(rep(0,lag),grow(x,Term-lag,growth))
    ans=zooreg(vals,start,end,frequency)
    return(ans)
  }
  x=c(rep(0,lag),x)
  if((length(x)!=Term)&warn) print('fullzoo warning -- number of values do not equal term,will pad zeros or truncate')
  if(length(x)<Term) {x=c(x,rep(0,Term-length(x)))}
  if(length(x)>Term) {x=x[1:Term]}
  ans=zooreg(x,start,end,frequency)
  return(ans)
}
#create an index to a zoo object with global variable defaults
fullzoo.ind=function(x=rep(0,Term),start=Start,end=End,frequency=Frequency) {
  if(length(x)!=Term) print('fullzoo warning -- number of values do not equal term')
  index(zooreg(x,start,end,frequency))
}
#create a yearmon for the fullzoo index
fullzoo.ym=function(x=rep(0,Term),start=Start,end=End,frequency=Frequency) {
  if(length(x)!=Term) print('fullzoo warning -- number of values do not equal term')
  as.yearmon(index(zooreg(x,start,end,frequency)))
}
#create zoo object with single value using global defaults for overall size
sparsezoo=function(v,date) {
  if(monthdiff(Start,date)<0) print('sparsezoo warning-date out of range')
  if((monthdiff(date,End)+1-length(v))<0) print('sparsezoo warning-date out of range')
  ans=fullzoo()
  v.z=zooreg(v,date,frequency=Frequency)
  si=match(index(v.z),fullzoo.ind())
  if(length(v)==1) {ans[si]=v}
  if(length(v)>1) {
    for (i in 1:(length(v))) {
      ans[si+i-1]=v[i]
    }    
  }
  ans
}
#assign in to a zoo object by month using global defaults
bymonzoo=function(v,month,growth=0) {
  ans=fullzoo()
  ans.ym=fullzoo.ym()
  if(is.numeric(month)) {
    month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")[month]
  }
  ans.ind=grep(month,ans.ym)
  if(length(ans.ind)==0) {
    print('bymonzoo warning -- no match, format for month is Jan, Feb, Mar, . . .')
    return(ans)
  }
  if(length(v)==1) {
    values=grow(v,length(ans.ind),growth)
    for (i in 1:length(ans.ind)) {
      ans[ans.ind[i]]=values[i]
    }
    return(ans)
  }
  if(length(v)>length(ans.ind)) print('bymonzoo warning -- too many values provided, not all used')
  if(length(v)<length(ans.ind)) print('bymonzoo warning -- not enough values, final value repeated')
  for (i in 1:(length(ans.ind))) {
    ans[ans.ind[i]]=v[min(i,length(v))]
  }
  return(ans)
}
# initiate empty data frames
initiate.cat=function(name) {
  ans=data.frame(fullzoo())
  colnames(ans)=name
  ans
}
# end a category by adding a total and removing the all zero 1st row
total.cat=function(x) {
  x$"Total"=rowSums(x)
  zoo(x,fullzoo.ind())
}
end.cat=function(x) {
  zoo(x,fullzoo.ind())
}
makeannualreport=function(x,...,orient="h",decim=0) {
  ans=round(toannual(merge(x,...)),decim)
  if (orient=="h") ans=t(ans)
  ans
}
allcases=function(cases) {
  ncol=length(cases)
  ans=matrix(0,prod(cases),ncol)
  ans[,ncol]=1:cases[ncol]
  if(ncol==1) return(ans)
  for (i in 1:(ncol-1)) {
    nprior=prod(cases[ncol:(ncol-i+1)])
    ans[,ncol-i]=rep(1:cases[ncol-i],rep(nprior,cases[ncol-i]))
  }
  return(ans)
  #old version below
  #ans=permutations(max(cases),length(cases),repeats.allowed=TRUE)
  #ans.ind=t(ans)<=cases
  #ans.ind=apply(ans.ind,2,all)
  #ans[ans.ind,]
}
casemat=function(x,...,names) {
  x=list(x,...)
  narg=length(x)
  chars=sapply(x,is.character)
  x.ind=as.matrix(allcases(sapply(x,length)))
  ncases=dim(x.ind)[1]
  ans.mat=matrix(" ",ncases,narg)
  ans2=list()
  for (i in 1:narg) {
    ans.mat[,i]=x[[i]][x.ind[,i]]
    ans2[[i]]=x[[i]][x.ind[,i]]
  }
  ans2=data.frame(ans2)
  colnames(ans2)=names
  for (i in 1:narg) {
    if (chars[i]) { 
    for (j in 1:ncases) {
       ans.mat[j,i]=paste("'",ans.mat[j,i],"'",sep="")
    }
  }
  }  
  for (i in 1:narg) {
    for (j in 1:ncases) {
      ans.mat[j,i]=paste(names[i],"=",ans.mat[j,i],sep="")
    }
  }
  ans=ans.mat[,1]
  if(narg==1) return(list(ans,ans2))
  for (i in 1:ncases) {
    for (j in 2:narg) {
      ans[i]=paste(ans[i],ans.mat[i,j],sep=",")      
    }
  }
  return(list(ans,ans2))
}
casemanager=function(fun,x,...,argnames,ansnames) {
  cases=casemat(x,...,names=argnames)
  case.str=cases[[1]]
  case.mat=cases[[2]]
  #for each case run fun and accumulate the answers in a list called answers
  answers=matrix(0,length(case.str),length(ansnames))
  for (i in 1:length(case.str)) {
  eval(parse(text=paste("temp=",fun,"(",case.str[i],")"))) 
  answers[i,]=as.vector(temp[[1]])
  }
  #organize the case data frame
  #data contract is that answer is a list and we take first element of the as the data.frame metrics
  answers=data.frame(answers)
  colnames(answers)=ansnames
  return(cbind(case.mat,answers))
}
#pe performance
#data is data frame zoo object, 1st column is cash flow (negative = draw), 2nd column is market index
pestats=function(cf,ind)  {
  ans=list()
    ans$tvpi=NA
    ans$irr=NA
    ans$pme=NA
    ans$pme.plus=NA
    ans$alpha=NA
    ans$ind.irr=NA
    ans$pme.wealth=NA
    ans$sb=NA
    ans$ss=NA
  if(length(cf)<=1)  return(ans)
  if(all(cf<=0)) return(ans)
  if(all(cf>=0)) return(ans)
  ans$tvpi=tvpi(cf)
  ans$irr=irr.z(cf,gips=TRUE)
  cf.neg=cf.pos=cf
  cf.pos[cf.pos<0]=0
  cf.neg[cf.neg>0]=0
  if(all(!(is.na(ind)))) {
  fvfactor=as.numeric(lastinvec(ind))/ind
  ans$pme=-sum(cf.pos*fvfactor)/sum(cf.neg*fvfactor)
  cf.fv=cf*fvfactor
  ans$alpha=log(1+irr.z(cf.fv,gips=TRUE))
  logpe.irr=log(1+ans$irr)
  logdm.irr=logpe.irr-ans$alpha
  ans$ind.irr=-1+exp(logdm.irr)
  ans$pme.wealthdiff=sum(cf.fv)
  sb=-cf.neg/ind
  if(length(which(cf.pos>0))>1) {
    sharesterm=as.numeric(lastinvec(cf.pos)/lastinvec(ind))
  stosell=sum(sb)-sharesterm
  cf.posxlast=cf.pos
  cf.posxlast[length(cf.pos)]=0
  scale=stosell/sum(cf.posxlast/ind)
  cf.pos.scaled=cf.pos*scale
  cf.pos.scaled[length(cf.pos)]=cf.pos[length(cf.pos)]
  ans$pme.plus=irr.z(mergesum.z(cf.pos.scaled,cf.neg),gips=TRUE) } else {
    cf.pos.scaled=cf.pos
    cf.pos.scaled[length(cf.pos)]=(sum(sb))*ind[length(ind)]
    ans$pme.plus=irr.z(mergesum.z(cf.pos.scaled,cf.neg),gips=TRUE)
  }
    
#   lndelt=-sum(fvfactor*cf)
#   cf.ln=cf
#   cf.ln[length(cf.ln)]=lndelt+cf.ln[length(cf.ln)]
#   ans$pme.ln=irr.z(cf.ln,gips=TRUE)
  ans$sb=sb
  ans$ss=cf.pos.scaled/ind
}
  
  return(ans)
}
pe.performance=function(cf.filename,bench.ticker.vec=c('^RUT','SPY'),
                        useallvals=FALSE,service="yahoo") {
  #function requires zoo and quantmod
  #cf.filename is the name of a file
  # in csv format, column one "Date" is date yyyy/mm/dd
  # column two is cash flow, labeled "Amount" , negative = draw
  # column three, labeled "Type", is a code "C" means a cash flow and "V" means a value, only the last value is used
  #bench.ticker.lst is a a vetor of ticker names available 
  # at yahoo finance
  #Returns a matrix with a column for each benchmark
  #Rows are 
  # IRR
  # TVPI
  # Realized %
  # PME (Kaplan)
  # Alpha (Griffiths)
  # PME+ (Capital Dynamics)
  if(is.data.frame(cf.filename)) {
  x=cf.filename} else {
  x=read.csv(file=cf.filename)
  }
# remove NAs
# if(any(is.na(x$Date))) {warning("Bad date formats");stop}
x$Date=as.Date(x$Date,format='%m/%d/%Y')
x$Amount=as.numeric(x$Amount)
badrows=which(is.na(x$Date)|is.na(x$Amount))
if (length(badrows>0)) {
  warning(paste("you have data causing NAs at row",paste(badrows,collapse=" ")))
  warning("removing these rows and running calculations")
  x=x[-badrows,]
}
x=subset(x,x$Amount!=0|x$Type=="V")


if(useallvals) {
  cf=aggregate(x$Amount,by=list(x$Date),sum)
  cf.wval=zoo(cf[,2],cf[,1])
  last=max(x$Date)
} else {
  x.cf=subset(x,x$Type=='C')
  x.v=subset(x,x$Type=='V')
  if(dim(x.v)[1]==0) {warning("No valuation");stop}
  tval=lastinvec(x.v$Amount)
  last=lastinvec(x.v$Date)
  tval=zoo(tval,last)
  #remove amounts after valuation date
  x.cf=subset(x.cf,x.cf$Date<=last)
  #aggregate by date and convert to zoo obect
  cf=aggregate(x.cf$Amount,by=list(x.cf$Date),sum)
  cf.noval=zoo(cf[,2],cf[,1])
  cf.wval=mergesum.z(cf.noval,tval)
}
#download benchmarks
first=min(x$Date)
bench.lst=list()
if (service == "yahoo") {
for (i in 1:length(bench.ticker.vec)) {
  rut=getSymbols("^RUT",src='yahoo',from=first-5,to=last+5,auto.assign=FALSE)
  bench.lst[[i]]=as.zoo(rut[,6])
}
} else {
  require(Rbbg)
  conn=blpConnect()
  first=as.Date(first,"%Y%m%d")
  last=as.Date(last,"%Y%m%d")
  for (i in 1:length(bench.ticker.vec)) {
    temp=bdh(conn,bench.ticker.vec[i],"PX_LAST",first-5,last+5)
    bench.lst[[i]]=zoo(temp$PX_LAST,as.Date(temp$date))
  }
  blpDisconnect(conn)
}
names(bench.lst)=bench.ticker.vec
#add missing days to benchmark
#fill NAs with neighboring values
days365=zooreg(rep(0,1+last-first),start=first,end=last)
for (i in 1:length(bench.ticker.vec)){
  b=merge(days365,bench.lst[[i]])
  b=b[,2]
  b=na.approx(b,na.rm=FALSE)
  bench.lst[[i]]=b
}
#match dates for benchmark reduction index
dind=match(index(cf.wval),index(bench.lst[[1]]))
#call pestats for each benchmark
ans=matrix(0,nrow=5,ncol=length(bench.lst))
for (i in 1:length(bench.lst)) {
  perf=pestats(cf.wval,bench.lst[[i]][dind])
  ans[,i]=c(perf$tvpi,100*perf$irr,perf$pme,100*perf$ind.irr,100*perf$alpha)
}
#build the return matrix
ans=round(ans,2)
colnames(ans)=bench.ticker.vec
rownames(ans)=c('TVPI','IRR','PME','Dollar Matched Index IRR','Direct Alpha')
return(ans)
}

pe.performance.return.index=function(cf.filename,bench.ticker.vec=c('^RUT','SPY'),
                        useallvals=FALSE,service="yahoo") {
  #function requires zoo and quantmod
  #cf.filename is the name of a file
  # in csv format, column one "Date" is date yyyy/mm/dd
  # column two is cash flow, labeled "Amount" , negative = draw
  # column three, labeled "Type", is a code "C" means a cash flow and "V" means a value, only the last value is used
  #bench.ticker.lst is a a vetor of ticker names available 
  # at yahoo finance
  #Returns a matrix with a column for each benchmark
  #Rows are 
  # IRR
  # TVPI
  # Realized %
  # PME (Kaplan)
  # Alpha (Griffiths)
  # PME+ (Capital Dynamics)
  if(is.data.frame(cf.filename)) {
    x=cf.filename} else {
      x=read.csv(file=cf.filename)
    }
  # remove NAs
  # if(any(is.na(x$Date))) {warning("Bad date formats");stop}
  x$Date=as.Date(x$Date,format='%m/%d/%Y')
  x$Amount=as.numeric(x$Amount)
  badrows=which(is.na(x$Date)|is.na(x$Amount))
  if (length(badrows>0)) {
    warning(paste("you have data causing NAs at row",paste(badrows,collapse=" ")))
    warning("removing these rows and running calculations")
    x=x[-badrows,]
  }
  x=subset(x,x$Amount!=0|x$Type=="V")
  
  
  if(useallvals) {
    cf=aggregate(x$Amount,by=list(x$Date),sum)
    cf.wval=zoo(cf[,2],cf[,1])
    last=max(x$Date)
  } else {
    x.cf=subset(x,x$Type=='C')
    x.v=subset(x,x$Type=='V')
    if(dim(x.v)[1]==0) {warning("No valuation");stop}
    tval=lastinvec(x.v$Amount)
    last=lastinvec(x.v$Date)
    tval=zoo(tval,last)
    #remove amounts after valuation date
    x.cf=subset(x.cf,x.cf$Date<=last)
    #aggregate by date and convert to zoo obect
    cf=aggregate(x.cf$Amount,by=list(x.cf$Date),sum)
    cf.noval=zoo(cf[,2],cf[,1])
    cf.wval=mergesum.z(cf.noval,tval)
  }
  #download benchmarks
  first=min(x$Date)
  bench.lst=list()
  if (service == "yahoo") {
    for (i in 1:length(bench.ticker.vec)) {
      rut=getSymbols("^RUT",src='yahoo',from=first-5,to=last+5,auto.assign=FALSE)
      bench.lst[[i]]=as.zoo(rut[,6])
    }
  } else {
    require(Rbbg)
    conn=blpConnect()
    first=as.Date(first,"%Y%m%d")
    last=as.Date(last,"%Y%m%d")
    for (i in 1:length(bench.ticker.vec)) {
      temp=bdh(conn,bench.ticker.vec[i],"PX_LAST",first-5,last+5)
      bench.lst[[i]]=zoo(temp$PX_LAST,as.Date(temp$date))
    }
    blpDisconnect(conn)
  }
  names(bench.lst)=bench.ticker.vec
  #add missing days to benchmark
  #fill NAs with neighboring values
  days365=zooreg(rep(0,1+last-first),start=first,end=last)
  for (i in 1:length(bench.ticker.vec)){
    b=merge(days365,bench.lst[[i]])
    b=b[,2]
    b=na.approx(b,na.rm=FALSE)
    bench.lst[[i]]=b
  }
  #match dates for benchmark reduction index
  dind=match(index(cf.wval),index(bench.lst[[1]]))
  #call pestats for each benchmark
  ans=matrix(0,nrow=5,ncol=length(bench.lst))
  for (i in 1:length(bench.lst)) {
    bench.lst[[i]]=bench.lst[[i]][dind]
    perf=pestats(cf.wval,bench.lst[[i]])
    ans[,i]=c(perf$tvpi,100*perf$irr,perf$pme,100*perf$ind.irr,100*perf$alpha)
  }
  #build the return matrix
  ans=round(ans,2)
  colnames(ans)=bench.ticker.vec
  rownames(ans)=c('TVPI','IRR','PME','Dollar Matched Index IRR','Direct Alpha')
  return(list(ans,bench.lst,cf.wval))
}

pe.performance.roll=function(cf,val,idx) {
#calculate pe performance statististics (TVPI, IRR, PME, PME+)
#cf is cash flow as a zoo object, length >= 1, must contain a negative value
#val is values as a zoo object, length >=1
#idx is index as a zoo object, must contain index value for all dates in cf and val
#returns a 5 column matrix
ans=matrix(0,nrow=length(val),col=6)
colnames(ans)=c('Date','IRR','TVPI','PME','PME+','Direct Alpha')
if(!is.zoo(cf)|!is.zoo(val)|!is.zoo(idx)) stop('pe.performance.roll arguments must all be zoo objects')
if (!all(index(cf) %in% index(idx))) stop('index object does not contain all dates in cash flow ')
if (!all(index(val) %in% index(idx))) stop('index object does not contain all dates in values')
for (i in length(val)) {
   v.i=val[i]
   date.i=index(v.i)
   cf.i=cf[index(cf)<=date.i]
   cf.i=mergesum.z(cf.i,v.i)
   dind=match(index(cf.i),index(idx))
   ans.i=pestats(cf.i,idx[dind])
   ans[i,1]=date.i
   ans[i,2]=ans.i$irr
   ans[i,3]=ans.i$tvpi
   ans[i,4]=ans.i$pme
   ans[i,5]=ans.i$pme.plus 
   ans[i,6]=ans.i$alpha
}
ans.l=list()
ans.l$irr=zoo(ans[,2],ans[,1])
ans.l$tvpi=zoo(ans[,3],ans[,1])
ans.l$pme=zoo(ans[,4],ans[,1])
ans.l$pme.plus=zoo(ans[,5],ans[,1])
return(ans.l)
}

#line of credit function
#args are beginning balance (default=0), interest rate and a cash flow
#calculates a line of credit assuming all negative cash is financed
#by the LOC. Interest is calculated on the average loan balance beginning 
#of period and end of period.  
#returns in a list the ending loan balance, interest expense and any leftover cash
#sign convention -- cash flows negative is cash out, positive cash in --
#for loan balance -- positive number reflects a loan balance
#uses solve to find answer to system of linear equations
lineofcredit=function(cf,ir,beg=0) {
  begbal=vector()
  begbal[1]=beg
  endbal=vector()
  intexp=vector()
  remain=vector()
  for (i in 1:length(cf)) {
    if(cf[i]>=((1+(ir/2))*begbal[i])) {
      intexp[i]=(ir/2)*begbal[i]
      endbal[i]=begbal[i+1]=0
      remain[i]=cf[i]-begbal[i]+endbal[i]-intexp[i]
    } else {
      coef=matrix(c(1,-1,-1,ir/2),byrow=TRUE,nrow=2)
      ys=c(cf[i]-begbal[i],-ir/2*begbal[i])
      ans=solve(coef,ys)
      intexp[i]=ans[1]
      endbal[i]=begbal[i+1]=ans[2]
      remain[i]=cf[i]-begbal[i]+endbal[i]-intexp[i]
    }
  }
  ans=list()
  ans$intexp=intexp
  ans$endbal=endbal
  ans$remain=remain
  return(ans)
}

## replacement for tseries get.hist.quote
get.hist.quote.fix=function (instrument = "^gdax", start, end, quote = c("Open", 
                                                      "High", "Low", "Close"), provider = c("yahoo", "oanda"), 
          method = NULL, origin = "1899-12-30", compression = "d", 
          retclass = c("zoo", "its", "ts"), quiet = FALSE, drop = FALSE) 
{
  if (missing(start)) 
    start <- "1991-01-02"
  if (missing(end)) 
    end <- format(Sys.Date() - 1, "%Y-%m-%d")
  provider <- match.arg(provider)
  retclass <- match.arg(retclass)
  start <- as.Date(start)
  end <- as.Date(end)
  if (is.null(method)) {
    method <- getOption("download.file.method")
    if (is.null(method)) 
      method <- "auto"
  }
  if (provider == "yahoo") {
    url <- paste("http://ichart.finance.yahoo.com/table.csv?s=", 
                 instrument, format(start, paste("&a=", as.character(as.numeric(format(start, 
                                                                                       "%m")) - 1), "&b=%d&c=%Y", sep = "")), format(end, 
                                                                                                                                     paste("&d=", as.character(as.numeric(format(end, 
                                                                                                                                                                                 "%m")) - 1), "&e=%d&f=%Y", sep = "")), "&g=", 
                 compression, "&q=q&y=0&z=", instrument, "&x=.csv", 
                 sep = "")
    destfile <- tempfile()
    i <- 1L
    repeat {
      status <- tryCatch(download.file(url, destfile, 
                                       method = method, quiet = quiet), error = identity)
      if (!inherits(status, "error") && (status == 0)) 
        break
      unlink(destfile)
      if (i >= 5L) {
        if (inherits(status, "error")) 
          stop(conditionMessage(status))
        else stop(sprintf("download error, status %d", 
                          status))
      }
      message("download error, retrying ...")
      i <- i + 1L
    }
    nlines <- length(count.fields(destfile, sep = "\n"))
    if (nlines == 1) {
      unlink(destfile)
      stop(paste("no data available for", instrument))
    }
    x <- read.table(destfile, header = TRUE, sep = ",", 
                    as.is = TRUE, fill = TRUE)
    x <- na.omit(x)
    if (nrow(x) >= 2L && x[1L, 1L] == x[2L, 1L]) {
      warning("first date duplicated, first instance omitted")
      x <- x[-1L, , drop = FALSE]
    }
    unlink(destfile)
    names(x) <- gsub("\\.", "", names(x))
    nser <- pmatch(quote, names(x)[-1]) + 1
    if (any(is.na(nser))) 
      stop("this quote is not available")
    n <- nrow(x)
    dat <- as.Date(as.character(x[, 1]), "%Y-%m-%d")
    if (!quiet && (dat[n] != start)) 
      cat(format(dat[n], "time series starts %Y-%m-%d\n"))
    if (!quiet && (dat[1] != end)) 
      cat(format(dat[1], "time series ends   %Y-%m-%d\n"))
    if (retclass == "ts") {
      jdat <- unclass(julian(dat, origin = as.Date(origin)))
      ind <- jdat - jdat[n] + 1
      y <- matrix(NA, nrow = max(ind), ncol = length(nser))
      y[ind, ] <- as.matrix(x[, nser, drop = FALSE])
      colnames(y) <- names(x)[nser]
      y <- y[, seq_along(nser), drop = drop]
      return(ts(y, start = jdat[n], end = jdat[1]))
    }
    else {
      x <- as.matrix(x[, nser, drop = FALSE])
      rownames(x) <- NULL
      y <- zoo(x, dat)
      y <- y[, seq_along(nser), drop = drop]
      if (retclass == "its") {
        if (inherits(tryCatch(getNamespace("its"), error = identity), 
                     "error")) 
          warning("package its could not be loaded: zoo series returned")
        else {
          index(y) <- as.POSIXct(index(y))
          y <- its::as.its(y)
        }
      }
      return(y)
    }
  }
  else if (provider == "oanda") {
    if (!missing(quote)) {
      warning("argument 'quote' ignored for provider 'oanda'")
    }
    if (!missing(compression)) {
      warning("argument 'compression' ignored for provider 'oanda'")
    }
    url <- paste("http://www.oanda.com/convert/fxhistory?lang=en&date1=", 
                 format(start, "%m"), "%2F", format(start, "%d"), 
                 "%2F", format(start, "%y"), "&date=", format(end, 
                                                              "%m"), "%2F", format(end, "%d"), "%2F", format(end, 
                                                                                                             "%y"), "&date_fmt=us&exch=", unlist(strsplit(instrument, 
                                                                                                                                                          split = "/"))[1], "&exch2=&expr=", unlist(strsplit(instrument, 
                                                                                                                                                                                                             split = "/"))[2], "&expr2=&margin_fixed=0&&SUBMIT=Get+Table&format=ASCII&redirected=1", 
                 sep = "")
    destfile <- tempfile()
    status <- download.file(url, destfile, method = method, 
                            quiet = quiet)
    if (status != 0) {
      unlink(destfile)
      stop(paste("download error, status", status))
    }
    x <- readLines(destfile, warn = quiet)
    unlink(destfile)
    if (length(grep("Sorry", x)) > 0) {
      msg <- unlist(strsplit(gsub("<[a-zA-Z0-9\\/]*>", 
                                  "", x[grep("Sorry", x)]), split = " "))
      msg <- paste(msg[msg != ""], collapse = " ")
      stop("Message from Oanda: ", msg)
    }
    first <- grep("<PRE>", x, fixed = TRUE)
    last <- grep("</PRE>", x, fixed = TRUE) - 1
    if ((length(first) != 1) || (length(last) != 1) || (last < 
                                                          first)) {
      stop(paste("no data available for", instrument))
    }
    x[first] <- sub(".*<PRE>", "", x[first])
    con <- textConnection(x[first:last])
    on.exit(close(con))
    x <- scan(con, what = list(character(), double()), quiet = TRUE)
    dat <- as.Date(x[[1]], format = "%m/%d/%Y")
    n <- length(dat)
    if (!quiet && (dat[1] != start)) 
      cat(format(dat[1], "time series starts %Y-%m-%d\n"))
    if (!quiet && (dat[n] != end)) 
      cat(format(dat[n], "time series ends   %Y-%m-%d\n"))
    if (retclass == "ts") {
      jdat <- unclass(julian(dat, origin = as.Date(origin)))
      ind <- jdat - jdat[1] + 1
      y <- rep.int(NA, max(ind))
      y[ind] <- x[[2]]
      return(ts(y, start = jdat[1], end = jdat[n]))
    }
    else {
      y <- zoo(x[[2]], dat)
      if (retclass == "its") {
        if (inherits(tryCatch(getNamespace("its"), error = identity), 
                     "error")) 
          warning("package its could not be loaded: zoo series returned")
        else {
          index(y) <- as.POSIXct(index(y))
          y <- its::as.its(y)
        }
      }
      return(y)
    }
  }
  else stop("provider not implemented")
}

blankfun=function(x,digits,addtxt,big.mark) {
  ifelse((x==0|is.na(x)|is.nan(x)),'   ',
         paste0(format(round(x,digits),big.mark=big.mark,nsmall=digits),addtxt)
         )
}
formatwblank=function(x,digits=0,addtxt=NULL,big.mark=',') {
  sapply(x,blankfun,digits,addtxt,big.mark)
}
grow.z=function(r,c,n,freq=1,start=1,retclass='zoo') {
  if('zoo'==retclass & start==1) {
    if(is.zoo(c)) start=time(c)
  }
  n=n*freq
  r=-1+(1+r)^(1/freq)
  ts=as.numeric(c)*(1+r)^(0:(n-1))
  if ('zoo'==retclass) {
    ts=zooreg(ts,start=start,freq=freq)
  }
  return(ts)
}

window.list=function(x,start,end) {
  ans=x
  for (i in 1:length(x)) {
    ans[[i]]=window(x[[i]],start=start,end=end)
  }
  return(ans)
}
addtotal.list=function(x) {
  total=do.call(zoosum,x)
  x$Total=total
  return(x)
}
merge0=function(...) {
  merge(...,fill=0)
}
make.table=function(x,by='year',time.horizontal=TRUE,fun=sum) {
  table=do.call(merge0,x)
  if (by=='year') table=aggregate(table,by=year(time(table)),fun)
  if (by=='quarter') table=aggregate(table,by=as.yearqtr(time(table)),fun)
  if (by=='month') table=aggregate(table,by=as.yearmon(time(table)),fun)
  timeline=time(table)
  table=as.matrix(table)
  rownames(table)=as.character(timeline)
  if(time.horizontal) table=t(table)
  return(table)
}
pfd.return=function(cf,int,freq=4,mdate=NA) {
  require(zoo)
  require(lubridate)
  #check validity of compounding period
  if(!freq%in%c(1,4,12)) stop('freq must be annual(1), quarterly(4) or monthly(12)')
  #check that cf is a zoo object
  if(!is.zoo(cf)) stop('cash flow must be zoo object')
  #ignore cash flow elements after the measurement date
  if(!is.na(mdate)) {
    cf=cf[time(cf)<=mdate]
    cf=zoosum(cf,zoo(0,mdate))
  } 
  tl=time(cf)
  if("Date"!=class(tl)) {
    tl=as.Date(tl)
    warning ('time index coerced to Date class')
  }
  if(freq==1) {
    time.ld=seq(year(tl[1]),year(lastinvec(tl)))
    time.ld=-1+as.Date(1+yearmon(time.ld))
    report.ld=.25+seq(as.yearqtr(tl[1]),as.yearqtr(lastinvec(tl)),by=.25)
    while(1!=quarter(lastinvec(report.ld))) {
      report.ld=c(report.ld,.25+lastinvec(report.ld))
    }
    report.ld=-1+as.Date(report.ld)
    if(is.na(mdate)) mdate=lastinvec(time.ld)
    if(lastinvec(time.ld)>mdate) {
      time.ld[length(time.ld)]=mdate
      report.ld=report.ld[report.ld<mdate]
      if(mdate!=lastinvec(report.ld)) report.ld=c(report.ld,mdate)
    }
  }
  if(freq==4) {
    time.ld=.25+seq(as.yearqtr(tl[1]),as.yearqtr(lastinvec(tl)),by=.25)
    time.ld=-1+as.Date(time.ld)
    if(is.na(mdate)) mdate=lastinvec(time.ld)
    if(lastinvec(time.ld)>mdate) {
      time.ld[length(time.ld)]=mdate
    }
    report.ld=time.ld
  }
  if(freq==12) {
    nmon=1+round(12*(as.yearmon(lastinvec(tl))-as.yearmon(tl[1])))
    time.ld=-1+as.Date(as.yearmon((tl[1])+months(1:nmon)))
    if(is.na(mdate)) mdate=lastinvec(time.ld)
    if(lastinvec(time.ld)>mdate) {
      time.ld[length(time.ld)]=mdate
    }
    report.ld=time.ld
  }
  acc.pref=current.pref=unreturned=zoo(0,report.ld)
  pref.paid=cap.paid=residual=0*cf
  time.comb=time(zoosum(cf,unreturned))
  unreturned.bal=current.bal=acc.bal=0
  n=length(time.comb)
  for (i in 1:n) {
    now=time.comb[i]
    ndays=0
    if(i!=1) ndays=as.numeric(time.comb[i]-time.comb[i-1])
    current.bal=current.bal+(int*ndays/365)*(acc.bal+unreturned.bal)
    if (now%in%tl) {
      cf.now=as.numeric(cf[now])
      if(cf.now>0) {
        split=waterfall(c(current.bal,acc.bal,unreturned.bal),cf.now)
        current.bal=current.bal-split[1]
        acc.bal=acc.bal-split[2]
        unreturned.bal=unreturned.bal-split[3]
        pref.paid[time.comb[i]]=split[1]+split[2]
        cap.paid[time.comb[i]]=split[3]
        residual[time.comb[i]]=split[4]
      } else {
        unreturned.bal=unreturned.bal-cf.now
      }
    }
    if(now%in%report.ld) {
      unreturned[now]=unreturned.bal
      acc.pref[now]=acc.bal
      current.pref[now]=current.bal
    }
    if(now%in%time.ld) {
      acc.bal=acc.bal+current.bal
      acc.pref[now]=acc.bal
      current.pref[now]=0
      current.bal=0
    }
  }
  ans=list()
  ans$current.pref=current.pref
  ans$acc.pref=acc.pref
  ans$unreturned=unreturned
  ans$pref.paid=pref.paid
  ans$cap.paid=cap.paid
  ans$residual=residual
  ans$cf=cf
  return(ans)
}

lsp <- function(package, all.names = FALSE, pattern) 
{
  package <- deparse(substitute(package))
  ls(
    pos = paste("package", package, sep = ":"), 
    all.names = all.names, 
    pattern = pattern
  )
}

coneplot=function(return,annual_expected_return,annual_standard_deviation,
                  periodicity=c('y','q','m'),ylabel=NULL) {
  #
  # return is a time series of monthly, quarterly or annual arithmetic return values
  # annual_expected_return is the annual expected return for the strategy
  # annual_standard_deviation is the expected annual volatility
  # periodicity is 'y', 'q' or 'm' indicating yearly, quarter or monthly data
  # ylabel is optional argument to provide a label for the y axis
  #
  # requires lubridate,zoo packages
  #
  n=length(return)
  months=0:n
  timelabel=months
  if (is.zoo(return)) {
    datehold=time(return)
    if(class(datehold)=="Date") {
      if (periodicity=='y') datehold=year(datehold)
      if (periodicity=='q') datehold=as.yearqtr(datehold)
      if (periodicity=='m') datehold=as.yearmon(datehold)
    }
    timelabel=as.numeric(datehold)
  }
  if(periodicity=='y') {
    expected_return=annual_expected_return
    standard_deviation=annual_standard_deviation
    if(is.zoo(return)) {
      timelabel=c(timelabel[1]-1,timelabel)
    }
    xlabel="Yearly Data"
  }
  if(periodicity=='q') {
    expected_return=-1+exp(log(1+annual_expected_return)/4)
    standard_deviation=annual_standard_deviation/sqrt(4)
    if(is.zoo(return)) {
      timelabel=c(timelabel[1]-(1/12),timelabel)
    }
    xlabel="Quarterly Data"
  }
  if(periodicity=='m') {
    expected_return=-1+exp(log(1+annual_expected_return)/12)
    standard_deviation=annual_standard_deviation/sqrt(12)
    if(is.zoo(return)) {
      timelabel=c(timelabel[1]-.25,timelabel)
    }
    xlabel="Monthly Data"
  }
  expected_gd=c(1,exp(cumsum(rep(log(1+expected_return),n))))
  actual_gd=c(1,exp(cumsum(log(1+coredata(return)))))
  one_sigma_up=expected_gd+(standard_deviation*months/sqrt(n))
  two_sigma_up=expected_gd+(2*standard_deviation*months/sqrt(n))
  one_sigma_down=expected_gd-(standard_deviation*months/sqrt(n))
  two_sigma_down=expected_gd-(2*standard_deviation*months/sqrt(n))
  conedf=data.frame(timelabel,expected_gd,actual_gd,one_sigma_up,two_sigma_up,
                    one_sigma_down,two_sigma_down)
  if (is.null(ylabel)) ylabel=names(x)[1]
  x=ggplot(conedf,(aes(x=timelabel)))+ylab(ylabel)+xlab(xlabel)+
    geom_line(linetype=2,colour='black',aes(y=expected_gd))+
    geom_line(colour='blue',aes(y=actual_gd))+
    geom_line(linetype=2,colour='orange',aes(y=one_sigma_up))+
    geom_line(linetype=2,colour='orange',aes(y=one_sigma_down))+
    geom_line(linetype=2,colour='red',aes(y=two_sigma_up))+
    geom_line(linetype=2,colour='red',aes(y=two_sigma_down))
  return(x)
}

unzoo=function(zoo.object,wide=TRUE,dataname="Value",datename="Date",category="Categories") {
  if(!is.zoo(zoo.object)) stop("not a zoo object")
  t=time(zoo.object)
  v=coredata(zoo.object)
  if(is.null(dim(zoo.object))) {
    ans=data.frame(t,v)
    colnames(ans)=c(datename,dataname)
    return(ans)
  }
  ans=as.data.frame(v)
  columnames=c(datename,colnames(v))
  ans=cbind(t,ans)
  colnames(ans)=columnames
  if(wide) return(ans)
  require(tidyr)
  colnames(ans)=c("Date",colnames(v))
  ans=gather(ans,cats,value,-Date)
  colnames(ans)=c(datename,category,dataname)
  return(ans)
}

incentfee.rev=function(nav,cf,cu,pref=.08,cry=.2,fre=4) {
  #
  # reverse engineer an incentive fee from following info
  # nav-- a one element time series of the investor net asset value
  # cf -- a time series of calls (negative) and distributions (positive), 
  #       values after date of nav are ignored
  # cu -- the catchup percent as a decimal
  # pref -- the preferred return rate as a decimal
  # cry -- the carry rate as a decimal
  # fre -- the periodicity of calculation, 4 for quarterly (the default and industry convention)
  #
  # returns a zoo data.frame object which is the incentive fee earned cumulatively through date of nav
  # only works for a one level structure with catchup
  # calculated on accrual and cash basis, both with european waterfall style calculations
  #
  # data validation
  if (!is.zoo(cf)) stop ("cf must be a zoo object")
  if (!is.zoo(nav)) stop ("nav must be a zoo object")
  if (cu<0|cu>1) stop ("cu must be 0<=cu<=1")
  if (pref<0|pref>1) stop ("pref must be 0<=pref<=1")
  if (cry<0|cry>1) stop ("cry must be 0<=cry<=1")
  if (cu!=0&cu<=cry) stop ("cu must be zero or > cry")
  if (!fre%in%c(1,4,12)) stop("freq must be 1, 4 or 12")
  cf=cf[time(cf)<=time(nav)]
  if(length(cf)==0) {
    warning ("no cash flows on or before nav date")
    return(NULL)
  }
  if(!any(cf<0)) {
    warning("no cash flows less than zero")
    return(NULL)
  }
  #calculate fee accrual basis
  #calculate preferred return dollars
  x=pfd.return(mergesum.z(cf,nav),int=pref,freq=fre,mdate=time(nav))
  prefearned=sum(x$pref.paid)+lastinvec(x$acc.pref)
  profitpaid=max(0,sum(x$pref.paid)+sum(x$residual)-lastinvec(x$unreturned))
  #calculate the incentive fee
  excess1=max(0,profitpaid-prefearned)
  #catch up layer
  lpsharecu=gpsharecu=0
  if(cu>0&excess1>0) {
    cuagg=prefearned*cry/(cu-cry)
    lpsharecu=(1-cu)*cuagg
    gpsharecu=cu*cuagg
    if(excess1<lpsharecu) {
      cumult=lpsharecu/excess1
      lpsharecu=lpsharecu*cumult
      gpsharecu=gpsharecu*cumult
    }
  }
  #parri passu layer
  excess2=excess1-lpsharecu
  incentfee=gpsharecu+(excess2*(-1+1/(1-cry)))
  
  #calculate fee paid w european waterfall
  x2=pfd.return(cf,int=pref,freq=fre,mdate=time(nav))
  prefearned=sum(x2$pref.paid)+lastinvec(x2$acc.pref)
  profitpaid=max(0,sum(x2$pref.paid)+sum(x2$residual)-lastinvec(x2$unreturned))
  #calculate the incentive fee
  excess1=max(0,profitpaid-prefearned)
  #catch up layer
  lpsharecu=gpsharecu=0
  if(cu>0&excess1>0) {
    cuagg=prefearned*cry/(cu-cry)
    lpsharecu=(1-cu)*cuagg
    gpsharecu=cu*cuagg
    if(excess1<lpsharecu) {
      cumult=lpsharecu/excess1
      lpsharecu=lpsharecu*cumult
      gpsharecu=gpsharecu*cumult
    }
  }
  #parri passu layer
  excess2=excess1-lpsharecu
  incentfee2=gpsharecu+(excess2*(-1+1/(1-cry)))
  
  #return answer
  return (data.frame(Accrual=zoo(incentfee,time(nav)),Cash=zoo(incentfee2,time(nav))))
}

fv.tr=function(CF,TR,ValDate=NULL) {  
  #
  # DESCRIPTION:
  # a function that calculates future value with time varying interest rates
  # 
  # ARGUMENTS:
  # CF -- a zoo object of cash flows
  # TR -- a zoo object which is a total return index reflecting the variable interest rate, 
  #       there must be an entry in TR for every date in CF and for the ValDate
  # ValDate -- the valuation date or dates, defaults to last date in cash flow if none provided,
  #             if earlier than last date in cash flows, then cash flows after this date are ignored,
  #             if earlier than first date in cash flows, returns NA
  #             defaults to the last date in the cash flow,
  #             must be in the same date format as CF and TR
  #
  # RESULT:
  # a zoo object of the same length as ValDate which is the value of the cash flows as if invested in the index
  # as of each valuation date
  #
  if(is.null(ValDate)) ValDate=tail(time(CF),1)
  if(!is.zoo(CF)) stop("CF not a zoo object")
  if(!is.zoo(TR)) stop("TR not a zoo object")
  if(any(!(ValDate%in%time(TR)))) stop("No TR for ValDate")
  if(any(!(time(CF)%in%time(TR)))) stop("No TR for CF")
  value=vector()
  for (i in 1:length(ValDate)) {
    if(ValDate[i]<min(time(CF))) value[i]=NA
    CF.i=CF[time(CF)<=ValDate[i]]
    TR.i=TR[time(TR)<=ValDate[i]]
    TR.ind=(as.numeric(tail(TR.i,1)))/TR.i
    value[i]=-sum(TR.ind[time(CF.i)]*CF.i)
  }
  return(zoo(value,ValDate))
}
IMD.palette=function(){
  return(c("#00619C","#DF8F26","#017769","#A61F2A",
           "#92D7EC","#FFCB78","#6DBCB6","#E27D84",
           "#886FA3","#4A2D7F","#A1A2A2","#F389B8"))
}
display.palette=function(pal=IMD.palette()) {
  pie(rep(1,length(pal)),labels=seq_along(pal),col=pal)
}
PME.palette=function(){
  return(c("#C05640", "#EDD170", "#1ECFD6", "#0878A4", "#003D73"))
}


#Copeland build from IMD\Karl\R projects\Functions.R

gd=function(x) {
  exp(cumsum(log(1+x)))
}

gdminus1=function(x,one=1) {
  temp=gd(x)
  temp-matrix(one,nrow=nrow(x),ncol=ncol(x),byrow=TRUE)
}

gdweight=function(x){
  temp=gd(x)
  temp/matrix(rowSums(temp),ncol=ncol(x),nrow=nrow(x))
}
gg=function(x,variable.name="Variable",value.name="Value") {
  melt(unzoo(x),id.vars="Date",
       variable.name=variable.name,
       value.name=value.name)
}
