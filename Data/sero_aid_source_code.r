############################################################
################ CREATING SEROAID OBJECT ###################
############################################################

create.seroaid.object<-function(age,seropos,group=NULL,age.min=1){

	##### Checking vector length #########
	
	cat('Creating SERO-AID data object\n')

	n.age<-length(age)
	n.seropos<-length(seropos)
	
	if(is.null(group)==T){
		
		group<-rep(1,n.age)
		
	}

	n.group<-length(group)

	if(n.age==n.group&n.seropos==n.group){
		
		cat('1. Data from',n.age,'individuals read;\n')

		##### Rounding age #########

		cat('2. Rounding age data to the closest integer;\n')	

		age<-round(age)

		##### Remove missing data #########
		
		keep.ids<-which(is.na(age)==F&is.na(seropos)==F&is.na(group)==F&age>=age.min)
	
		cat('3. Keeping data from ',length(keep.ids),' individuals with complete information and age equal to or greater than ',age.min,';\n',sep='')			
		age.filtered<-age[keep.ids]
		
		seropos.filtered<-seropos[keep.ids]
		
		group.filtered<-group[keep.ids]
		
		out<-list(age=age.filtered,seropos=seropos.filtered,group=as.character(group.filtered),class='sero-aid')
		
		cat('4. Data object successfully saved.\n')
		
	} else {
		
		cat("ERROR: Vectors 'age', 'seropos' or 'group' have different lengths\n")
		out<-NA
	}	
	
	return(out)
	
}

age.profile<-function(seroaid.object,analysis='overall',lag=0.05,conf.level=0.95){
	
	if(seroaid.object$class=='sero-aid'){
	
		age<-seroaid.object$age
		sero.pos<-seroaid.object$seropos
		groups<-sort(unique(seroaid.object$group))
		
		if(analysis=='overall'){
					
			output<-age.profile.group(age,sero.pos,lag,conf.level)
			group<-rep('overall',dim(output)[1])
			
			output<-data.frame(group=group,output)
			output<-list(age.profiles=output,analysis='overall',class='sero-aid')
		
		} else {	

			if(analysis=='split'){
				
				output<-c()
				group<-c()
				
				for(i in groups){

					lista<-which(seroaid.object$group==i)	
			
					results<-age.profile.group(age[lista],sero.pos[lista],lag, conf.level)
						
					group<-c(group,rep(i,dim(results)[1]))
			
					output<-rbind(output,results)

				}
				output<-data.frame(group=group,output)
				output<-list(age.profiles=output,analysis='split',class='sero-aid')
							
			} else { 	

				if(analysis!='overall'&analysis!='split'){
		
					exist.group<-sum(groups==analysis)
		
					if(exist.group==1){
	
						lista<-which(seroaid.object$group==analysis)	
						output<-age.profile.group(age[lista],sero.pos[lista],lag, conf.level)
						group<-rep(analysis,dim(output)[1])
						output<-data.frame(group=group,output)
		
						output<-list(age.profiles=output,analysis=analysis,class='sero-aid')
							
					} else {
			
						cat("ERROR: Group name does not exist in the data set!\t")
					
					}
				}	
			}		
		}			
	} else {
		
		cat("ERROR: Data object does not belong to sero-aid class. Please use 'create.seroaid.object' and create your data object first.")
		output<-NA
	
	}
	return(output)
	
}

##########################################
################ PLOTS ###################
##########################################

age.profile.group<-function(age,sero.pos,lag,conf.level){

	list.prop<-seq(lag,1,lag)

	list.quantiles<-unique(round(quantile(age,list.prop)))

	num.points<-length(list.quantiles)
		
	output<-matrix(NA,ncol=7,nrow=num.points)	

	list.aux<-which(age<=list.quantiles[1])
	output[1,1]<-min(age[list.aux])
	output[1,2]<-max(age[list.aux])
	output[1,3]<-length(sero.pos[list.aux])
	output[1,4]<-sum(sero.pos[list.aux])	
	output[1,5]<-mean(sero.pos[list.aux])
	output[1,6:7]<-binom.test(sum(sero.pos[list.aux]),length(list.aux),conf.level=conf.level)$conf.int
	
		for(i in 2:num.points){
	
			list.aux<-which(age>list.quantiles[i-1]&age<=list.quantiles[i])
			output[i,1]<-min(age[list.aux])
			output[i,2]<-max(age[list.aux])
			output[i,3]<-length(sero.pos[list.aux])
			output[i,4]<-sum(sero.pos[list.aux])	
			output[i,5]<-mean(sero.pos[list.aux])
			output[i,6:7]<-binom.test(sum(sero.pos[list.aux]),length(list.aux),conf.level=conf.level)$conf.int
	
		}
		
		output<-data.frame(output)
		colnames(output)<-c('age.lower','age.upper','n','n.seropos','sero.prev','lower','upper')
		return(output)

}


plot.seroaid<-function(sero.obs,sero.model=NULL,xlab='age (in years)',ylab='seroprevalence',ylim=NULL,xlim=NULL,lwd=NULL,lege=NULL,col.obs=NULL,col.exp=NULL,main=NULL){

	#### Defining defaults ######	

	if(is.null(xlim)==T)xlim<-c(0,10*ceiling(max((sero.obs$age.profiles)$age.upper/10)))
	if(is.null(ylim)==T)ylim<-c(0,1)	
	if(is.null(lwd)==T)lwd<-2	
	if(is.null(main)==T)main<-sero.obs$analysis
	if(is.null(col.obs)==T){
		
		num.col<-length(unique((sero.obs$age.profiles)$group))
		col.obs<-rainbow(num.col)
	}	
	if(is.null(col.exp)==T){
		
		num.col<-length(unique((sero.obs$age.profiles)$group))
		col.exp<-rainbow(num.col)
	}
		
	#### Plotting ###############	
		
	if(is.null(sero.obs)==F){

		if(sero.obs$class=='sero-aid'){
	
			plot(sero.obs$age.lower,sero.obs$seroprev,type='n',axes=F,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
			
			box()
		
			axis(1,seq(xlim[1],xlim[2],1),rep('',1+(xlim[2]-xlim[1])),las=1,lwd.ticks=0.5,tcl=-0.25)		
			axis(1,seq(xlim[1],xlim[2],10),las=1,lwd.ticks=1)
			axis(2,seq(ylim[1],ylim[2],0.01),rep('',1+(ylim[2]-ylim[1])/0.01),las=1,lwd.ticks=0.5,tcl=-0.25)
			axis(2,seq(ylim[1],ylim[2],0.1),las=1,lwd.ticks=1)
					
			if(is.null(lege)==F){
				
				groups<-colnames(sero.model$expected.seroprevalence)[-1]
				num.groups<-length(groups)
				colours1<-c()
				colours2<-c()
				groups1<-c()
				pch<-c()
				lwd2<-c()
				for(i in 1:num.groups){

					colours1<-c(colours1,col.obs[i],NA)
					colours2<-c(colours2,NA,col.exp[i])
					groups1<-c(groups1,paste(groups[i],'obs',sep=' - '),paste(groups[i],'exp',sep=' - '))
					lwd2<-c(lwd2,NA,1)
					pch<-c(pch,21,NA)
				}
				
				legend(lege[1],lege[2],groups1,lwd=lwd,pch=pch,pt.bg=colours1,col=colours2)
				
			}
			
		
			if(is.null(sero.model)==F){

				k<-1			
				for(i in 2:dim(sero.model$expected.seroprevalence)[2]){		
				
					lines(c(0,sero.model$expected.seroprevalence[,1]),c(0,sero.model$expected.seroprevalence[,i]),lwd=lwd,col=col.exp[k])
				
					k<-k+1
								
				}		
			
			}
			
			k<-1

			for(i in sort(unique((sero.obs$age.profiles)$group))){		
				
				aux<-which((sero.obs$age.profiles)$group==i)	
				
				#lines(c(0,(sero.obs$age.profiles)$age.upper[aux]),c(0,(sero.obs$age.profiles)$sero.prev[aux]),type='s',lwd=lwd,col=col.obs[k])
				points((sero.obs$age.profiles)$age.upper[aux],(sero.obs$age.profiles)$sero.prev[aux],pch=19,col=col.obs[k])
				k<-k+1
								
			}
		 
		 } else { 
		 	
		 	cat("ERROR: Data object does not belong to sero-aid class. Please use 'create.seroaid.object' and create your data object first.\n")
		 	
		 }	
	} else cat('ERROR: please create observed age profile\n')	

}

####################################################################################
################ ESTIMATION OF SIMPLE REVERSIBLE CATALYTIC MODEL ###################
####################################################################################

rcm.constant.analysis<-function(seroaid.object,analysis='overall',n.boot=NULL,conf.level=0.95,num.ini=1){
	
	#### errors #####
	
	error.class<-seroaid.object$class!='sero-aid'
	
	if(error.class==0){
		
		groups<-unique(seroaid.object$group)
		analyses<-c(groups,'overall','split','shared-scr-srr')
		error.analysis<-(1-sum(analysis==analysis))	
		
		###### estimation of the models #######
		
		if(error.analysis==0){
	
			if(analysis=='split'){
				
				output<-rcm.constant.single.group(seroaid.object$age,seroaid.object$seropos,seroaid.object$group,analysis,conf.level=conf.level)
				
			}
			if(analysis=='overall'){
				
				seroaid.object$group<-rep('overall',length(seroaid.object$group))
				output<-rcm.constant.single.group(seroaid.object$age,seroaid.object$seropos,seroaid.object$group,analysis,conf.level=conf.level)
				
			}

			if(analysis=='shared-scr-srr'){
				
				output<-rcm.constant.estimation.shared.scr.srr(seroaid.object,analysis,conf.level=conf.level)
								
			}
	
			if(analysis=='shared-scr'){
					
				output<-rcm.constant.estimation.shared.scr(seroaid.object,analysis,n.boot,conf.level=conf.level,num.ini=num.ini)
							
			}
	
			if(analysis=='shared-srr'){
				
				output<-rcm.constant.estimation.shared.srr(seroaid.object,analysis,n.boot,conf.level=conf.level,num.ini=num.ini)
								
			}
	
	
			if(sum(analysis==groups)==1){
				
				aux<-which(seroaid.object$group==analysis)
				output<-rcm.constant.single.group(seroaid.object$age[aux],seroaid.object$seropos[aux],seroaid.object$group[aux],analysis,conf.level=conf.level)
			}
			return(output)
										
		} else {
			
			cat('ERROR: Unknown option in "analysis" option\n')	
		}	
		
	} else {
		
		cat("ERROR: Data object does not belong to sero-aid class. Please use 'create.seroaid.object' and create your data object first.")
		
	}
}


####################################################
################ SIMPLE ANALYSIS ###################
####################################################

seromodel.constant<-function(x,log.lambda,log.rho)exp(log.lambda)/(exp(log.lambda)+exp(log.rho))*(1-exp(-(exp(log.lambda)+exp(log.rho))*x))

loglik.rcm.constant.single.group<-function(par,age,seropos){
	
	tabela.aux<-table(age,factor(seropos,levels=c(0,1)))
		
	n1<-tabela.aux[,2]
	n0<-tabela.aux[,1]
	t<-as.numeric(rownames(tabela.aux))
	
	p.binom<-as.numeric(sapply(t,seromodel.constant,log.lambda=par[1],log.rho=par[2]))
	data<-cbind(n1,n1+n0)
	data<-cbind(data,p.binom)
	loglik<-as.numeric(apply(data,1,function(x)dbinom(x=x[1],size=x[2],prob=x[3],log=T)))
	
	loglik.total<-sum(loglik)
	return(loglik.total)
}



rcm.constant.single.group<-function(age,seropos,group,analysis,par=NULL,conf.level){
		
	output<-c()
	
	estimates<-c()
		
	age.max<-max(age)
	
	fitted.values<-c(1:age.max)		

	loglik<-0

	for(g in sort(unique(group))){
		
		aux<-which(group==g)
	
		if(is.null(par)==T)par<-runif(2,-3,-1)
	
		sol<-optim(par,loglik.rcm.constant.single.group,age=age[aux],seropos=seropos[aux],control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	
		lambda1<-exp(sol$par[1])
		rho1<-exp(sol$par[2])

		prob.exp<-sapply(1:age.max,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda1,rho=rho1)	

		loglik<-loglik+sol$value

		fitted.values<-cbind(fitted.values,prob.exp)
		
		conf.int.rho<-rcm.constant.single.group.conf.int.rho(age[aux],seropos[aux],log(rho1),sol$value,conf.level)$rho
				
		conf.int.lambda<-rcm.constant.single.group.conf.int.lambda(age[aux],seropos[aux],log(lambda1),sol$value,conf.level)$lambda
		
		estimates<-rbind(estimates,c(lambda1,conf.int.lambda[1],conf.int.lambda[2],rho1,conf.int.rho[1],conf.int.rho[2]))

	}
	
	estimates<-data.frame(group=sort(unique(group)),round(estimates,4))

	colnames(estimates)<-c('group','lambda.est','lambda.lower','lambda.upper','rho.est','rho.lower','rho.upper')
		
	colnames(fitted.values)<-c('age',sort(unique(group)))
	
	output<-list(loglik.total=loglik,estimates=estimates,expected.seroprevalence=fitted.values,df=dim(estimates)[1]*2,model='M1',analysis=analysis)
		
	return(output)
	
}

####################################################
############ CONFIDENCE INTERVAL ###################
####################################################

rcm.constant.single.group.estimating.lambda.fixed.rho<-function(log.rho,age,seropos){
	
	f1<-function(log.lambda,log.rho,age,seropos)loglik.rcm.constant.single.group(par=c(log.lambda,log.rho),age,seropos)
	
	sol<-optimize(f1,c(-15,10),log.rho=log.rho,age=age,seropos=seropos,maximum=T)
	
	return(sol)
	
}	

rcm.constant.single.group.estimating.rho.fixed.lambda<-function(log.lambda,age,seropos){
	
	f1<-function(log.lambda,log.rho,age,seropos)loglik.rcm.constant.single.group(par=c(log.lambda,log.rho),age,seropos)
	
	sol<-optimize(f1,c(-15,10),log.lambda=log.lambda,age=age,seropos=seropos,maximum=T)
		
	return(sol)
	
}


	
rcm.constant.single.group.conf.int.rho<-function(age,seropos,log.rho,loglik.1,conf.level=0.95){

	f2<-function(log.rho,age,seropos,loglik.0)rcm.constant.single.group.estimating.lambda.fixed.rho(log.rho,age,seropos)$objective-loglik.0

	p<-qchisq(conf.level^(1/2),1)
	
	loglik.0<-p/(-2)+loglik.1
		
	##### Calculating lower bound ######
			
	sol<-rcm.constant.single.group.estimating.lambda.fixed.rho(-25,age,seropos)
			
	if(sol$objective>loglik.0){
		
		lower.rho<-0
				
	} else {
						
		sol.rho<-uniroot(f2,c(-25,log.rho),age=age,seropos=seropos,loglik.0=loglik.0)	
		lower.rho<-exp(sol.rho$root)
	
	}
	
	##### Calculating Upper bound ######
		
	sol<-rcm.constant.single.group.estimating.lambda.fixed.rho(10,age,seropos)
			
	if(sol$objective>loglik.0){
		
		upper.rho<-Inf
		
	} else {

		sol.rho<-uniroot(f2,c(log.rho,10),age=age,seropos=seropos,loglik.0=loglik.0)	
		upper.rho<-exp(sol.rho$root)

	}	
	return(list(rho=c(lower.rho,upper.rho)))
	
}


rcm.constant.single.group.conf.int.lambda<-function(age,seropos,log.lambda,loglik.1,conf.level=0.95){

	f2<-function(log.lambda,age,seropos,loglik.0)rcm.constant.single.group.estimating.rho.fixed.lambda(log.lambda,age,seropos)$objective-loglik.0

	p<-qchisq(conf.level^(1/2),1)
	
	loglik.0<-p/(-2)+loglik.1
		
	##### Calculating lower bound ######
			
	sol<-rcm.constant.single.group.estimating.rho.fixed.lambda(-25,age,seropos)
		
	if(sol$objective>loglik.0){
		
		lower.lambda<-0
				
	} else {
						
		sol.lambda<-uniroot(f2,c(-25,log.lambda),age=age,seropos=seropos,loglik.0=loglik.0)	
		lower.lambda<-exp(sol.lambda$root)
				
	}
	
	##### Calculating Upper bound ######
		
	sol<-rcm.constant.single.group.estimating.rho.fixed.lambda(10,age,seropos)
		
	if(sol$objective>loglik.0){
		
		upper.lambda<-Inf
		
	} else {

		sol.lambda<-uniroot(f2,c(log.lambda,10),age=age,seropos=seropos,loglik.0=loglik.0)	
		upper.lambda<-exp(sol.lambda$root)

	}	
	return(list(lambda=c(lower.lambda,upper.lambda)))
	
}

###############################################
################ SHARED SRR/SCR ###############
###############################################

rcm.constant.estimation.shared.scr.srr<-function(seroaid.object,analysis,conf.level){

	cat('Generating initial conditions\n')
	
	output<-rcm.constant.single.group(seroaid.object$age,seroaid.object$seropos,seroaid.object$group,'split',conf.level=conf.level)
	
	estimates.lambda<-log(runif(1,0,0.5))
	estimates.rho<-log(runif(1,0,0.5))

	par=c(estimates.lambda,estimates.rho)

	cat('Performing MLE\n')

	groups<-sort(unique(seroaid.object$group))
	num.group<-length(groups)
	
	seroaid.object2<-seroaid.object
	seroaid.object2$group<-rep('overall',length(seroaid.object2$group))
	
	output<-rcm.constant.single.group(seroaid.object2$age,seroaid.object2$seropos,seroaid.object2$group,analysis,par=par,conf.level=conf.level)
	
	estimates<-output$estimates
	seroprev<-output$expected.seroprevalence[,2]
	
	#if(estimates[1,5]==0)estimates[1,5]<-exp(-100)

	par<-c(log(estimates[1,2]),log(estimates[1,5]))
		
	output$loglik.total<-0
	output$estimates<-c()
	output$expected.seroprevalence<-output$expected.seroprevalence[,1]

	for(g in 1:num.group){
		
		aux<-which(seroaid.object$group==groups[g])
		
		output$loglik.total<-output$loglik.total+loglik.rcm.constant.single.group(par=par,seroaid.object2$age[aux],seroaid.object2$seropos[aux])	
		output$estimates<-rbind(output$estimates,estimates)
		output$expected.seroprevalence<-cbind(output$expected.seroprevalence,seroprev)
	
	}
	
	output$estimates[,1]<-groups

	colnames(output$expected.seroprevalence)<-c('age',groups)
		
	return(output)
}

###############################################
################ SHARED SRR ###################
###############################################

loglik.rcm.constant.shared.srr<-function(seroaid.object,par){
	
	loglik<-0
	groups<-unique(seroaid.object$group)
	num.groups<-length(groups)
	age.max<-max(seroaid.object$age)

	log.rho1<-par[num.groups+1]

	for(i in 1:num.groups){
		
		aux<-which(seroaid.object$group==groups[i])
		
		log.lambda1<-par[i]
		
		prob.exp<-sapply(1:age.max,function(x,log.lambda,log.rho)exp(log.lambda)/(exp(log.lambda)+exp(log.rho))*(1-exp(-(exp(log.lambda)+exp(log.rho))*x)),log.lambda=log.lambda1,log.rho=log.rho1)
		
		tabela<-table(factor(seroaid.object$age[aux],levels=1:age.max),factor(seroaid.object$seropos[aux],levels=c(0,1)))
		
		#print(tabela)	
				
		loglik<-loglik+sum(apply(cbind(tabela[,2],rowSums(tabela),prob.exp),1,function(x)dbinom(x=x[1],size=x[2],prob=x[3],log=T)))
				
	}
	return(loglik)
}

rcm.constant.estimation.shared.srr<-function(seroaid.object,analysis,n.boot=NULL,conf.level,num.ini){
	
	num.groups<-length(unique(seroaid.object$group))
	groups<-unique(seroaid.object$group)
	aux.order<-order(unique(seroaid.object$group))
	age.max<-max(seroaid.object$age)
	
	#print(groups[aux.order])
	
	#### Generating initial conditions #### 
	
	#cat('Generating initial conditions\n')
	
	output<-rcm.constant.single.group(seroaid.object$age,seroaid.object$seropos,seroaid.object$group,'split',conf.level=conf.level)
		
	#print(output$estimates)
	#cat('Performing MLE\n\n')
	
	loglik.max<-(-10^6)
	
	for(n in 1:num.ini){
	
		cat('Using initial value',n,'out of',num.ini,'\n')	
		estimates.lambda<-log(output$estimates[aux.order,2])+runif(num.groups,-0.1,0.1)
		estimates.rho<-log(mean(output$estimates[aux.order,5]))
		par<-c(estimates.lambda,estimates.rho)
	
		sol<-optim(par,loglik.rcm.constant.shared.srr,seroaid.object=seroaid.object,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	
		if(sol$value>loglik.max){
			
			loglik.max<-sol$value
			sol.best<-sol
		}
	}
		
	sol<-sol.best
		
	fitted.values<-1:age.max	
	estimates<-c()

	for(i in 1:num.groups){
		
		lambda1<-exp(sol$par[i])
		rho1<-exp(sol$par[num.groups+1])

		estimates<-rbind(estimates,c(lambda1,NA,NA,rho1,NA,NA))
		prob.exp<-sapply(1:age.max,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda1,rho=rho1)	

		fitted.values<-cbind(fitted.values,prob.exp)
	}
	
	if(is.null(n.boot)==T){
		
		#cat('Calculating CI by inverting LRT\n')

		#conf.int.lambda<-rcm.constant.estimation.shared.srr.conf.int.lambda(seroaid.object,sol$par[1:num.groups],sol$value,conf.level)

		#conf.int.rho<-rcm.constant.estimation.shared.srr.conf.int.rho(seroaid.object,sol$par[num.groups+1],sol$value,estimates.lambda,conf.level)
				
		#estimates[,2:3]<-conf.int.lambda$lambda
		#estimates[,5:6]<-conf.int.rho$rho
	
	} else {
		
		#cat('Calculating CI by bootstrap\n')
				
		estimates<-rcm.constant.estimation.shared.srr.conf.int.boot(seroaid.object,estimates,n.boot,conf.level)
		

	}
		
	estimates<-data.frame(group=groups,round(estimates,4))

	colnames(estimates)<-c('group','lambda.est','lambda.lower','lambda.upper','rho.est','rho.lower','rho.upper')
	
	colnames(fitted.values)<-c('age',groups)
	
	fitted.values<-cbind(fitted.values[,1],(fitted.values[,2:(num.groups+1)])[,aux.order])
		
	output<-list(loglik.total=sol$value,estimates=estimates[aux.order,],expected.seroprevalence= fitted.values,df=dim(estimates)[1]+1,model='M1',analysis=analysis)
	
}

###################################################################
################## Confidence intervals (LRT) #####################
###################################################################

############################# 1. Rho ###############################

loglik.rcm.constant.shared.srr.fixed.rho<-function(seroaid.object,log.lambda,log.rho)loglik.rcm.constant.shared.srr(seroaid.object,par=c(log.lambda,log.rho))		

rcm.constant.shared.srr.estimating.lambda.fixed.rho<-function(seroaid.object,log.rho,lambda.ini){

	sol<-optim(par=lambda.ini,loglik.rcm.constant.shared.srr.fixed.rho,seroaid.object=seroaid.object,log.rho=log.rho,control=list(fnscale=-1,maxit=100000,pgtol=1E-10))
	return(sol)
	
}

rcm.constant.estimation.shared.srr.conf.int.rho<-function(seroaid.object,log.rho,loglik.1,lambda.est,conf.level){

	#cat('Shared Rho\n',sep='')

	f2<-function(log.rho,seroaid.object,loglik.0,lambda.ini)rcm.constant.shared.srr.estimating.lambda.fixed.rho(seroaid.object,log.rho,lambda.ini)$value-loglik.0
	
	num.par<-length(unique(seroaid.object$group))+1
	
	p<-qchisq(conf.level^(1/num.par),1)
				
	loglik.0<-p/(-2)+loglik.1
		
	##### Calculating lower bound ######
	
	sol<-rcm.constant.shared.srr.estimating.lambda.fixed.rho(seroaid.object,-10,lambda.est)
		
	if(sol$value>loglik.0|log.rho<(-10)){
		
		lower.rho<-0
				
	} else {
						
		sol.rho<-uniroot(f2,c(-10,log.rho),seroaid.object=seroaid.object,loglik.0=loglik.0,lambda.ini=lambda.est)
		lower.rho<-exp(sol.rho$root)

	}
		
	##### Calculating Upper bound ######
		
	sol<-rcm.constant.shared.srr.estimating.lambda.fixed.rho(seroaid.object,10,lambda.est)
		
	if(sol$value>loglik.0|log.rho>10){		
		
		upper.rho<-Inf
		
	} else {

		sol.rho<-uniroot(f2,c(log.rho,10),seroaid.object=seroaid.object,loglik.0=loglik.0,lambda.ini=lambda.est)
		upper.rho<-exp(sol.rho$root)
	}	
	
	rho<-c()
	for(i in 1:(num.par-1))rho<-rbind(rho,c(lower.rho,upper.rho))
	
	return(list(rho=rho))
	
}

########################### 2. Lambda ###################################

loglik.rcm.constant.shared.srr.fixed.lambda<-function(seroaid.object,param,log.lambda1,j){
	
	g<-length(unique(seroaid.object$group))

	if(j==1)all.param<-c(log.lambda1,param[1:g])
	
	if(j==g)all.param<-c(param[1:(g-1)],log.lambda1,param[g])
	
	if(j!=1&j!=g)all.param<-c(param[1:(j-1)],log.lambda1,param[j:g])
		
	return(loglik.rcm.constant.shared.srr(seroaid.object,par=all.param))		

}


rcm.constant.shared.srr.estimating.rho.fixed.lambda<-function(seroaid.object,log.lambda1,j){
	
	g<-length(unique(seroaid.object$group))	
	sol<-optim(par=rep(-2,g), loglik.rcm.constant.shared.srr.fixed.lambda,seroaid.object=seroaid.object, log.lambda1=log.lambda1,j=j,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	return(sol)
	
}


rcm.constant.estimation.shared.srr.conf.int.lambda<-function(seroaid.object,log.lambda,loglik.1,conf.level){

	f2<-function(log.lambda1,seroaid.object,loglik.0,j)rcm.constant.shared.srr.estimating.rho.fixed.lambda(seroaid.object,log.lambda1,j)$value-loglik.0

	g<-length(unique(seroaid.object$group))

	num.par<-g+1
		
	p<-qchisq(conf.level^(1/num.par),1)
				
	loglik.0<-p/(-2)+loglik.1

	lambda<-c()

	for(i in 1:g){
		
		#cat('Lambda_',i,'\n',sep='')
		
		estimate<-log.lambda[i]
			
		##### Calculating lower bound ######

		sol<-rcm.constant.shared.srr.estimating.rho.fixed.lambda(seroaid.object,-10,i)
		
		#print(loglik.0)
		#print(sol)
			
		if(sol$value>loglik.0){
		
			lower.lambda<-0
				
		} else {
						
			sol.lambda<-uniroot(f2,c(-10,estimate),seroaid.object=seroaid.object,loglik.0=loglik.0,j=i)	
			lower.lambda<-exp(sol.lambda$root)
	
		}
	

		##### Calculating Upper bound ######
		
		sol<-rcm.constant.shared.srr.estimating.rho.fixed.lambda(seroaid.object,10,i)
				
		if(sol$value>loglik.0){
		
			upper.lambda<-(+Inf)
		
		} else {

			sol.lambda<-uniroot(f2,c(estimate,10),seroaid.object=seroaid.object,loglik.0=loglik.0,j=i)	
			upper.lambda<-exp(sol.lambda$root)
		}	
	
		lambda<-rbind(lambda,c(lower.lambda,upper.lambda))
	}
	
	return(list(lambda=lambda))
	
}

###################################################################
################## Confidence intervals (BOOTSTRAP) ###############
###################################################################

rcm.constant.estimation.shared.srr.conf.int.boot<-function(seroaid.object,estimates,n.boot,conf.level){
	
	num.groups<-length(unique(seroaid.object$group))
	
	estimates.bootstrap<-c()
	
	for(i in 1:n.boot){
		
		cat('Analysing bootstrap sample',i,'\n')
		
		seroaid.object.boot<-seroaid.object
			
		groups<-unique(seroaid.object$group)
		
		for(g in 1:num.groups){
		
			aux<-which(seroaid.object.boot$group==groups[g])
			
			lambda<-estimates[g,1]
			rho<-estimates[g,4]
			
			seroaid.object.boot$seropos[aux]<-sapply(seroaid.object.boot$age[aux],function(t,lambda,rho)rbinom(1,1,prob=(lambda/(lambda+rho)*(1-exp(-(lambda+rho)*t)))),lambda=lambda,rho=rho)

			
		}	
		
		#### Generating initial conditions #### 

		par<-log(c(estimates[,1],estimates[1,4]))
			
		sol<-optim(par,loglik.rcm.constant.shared.srr,seroaid.object=seroaid.object.boot,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
		
		estimates.bootstrap<-rbind(estimates.bootstrap,sol$par)
		
	}
	
	for(g in 1:num.groups){
		
		estimates[g,2]<-exp(quantile(estimates.bootstrap[,g],(1-conf.level)/2))
		estimates[g,3]<-exp(quantile(estimates.bootstrap[,g],(conf.level+1)/2))
		estimates[g,5]<-exp(quantile(estimates.bootstrap[,num.groups+1],(1-conf.level)/2))
		estimates[g,6]<-exp(quantile(estimates.bootstrap[,num.groups+1],(conf.level+1)/2))
	
	}
	
	return(estimates)
	
}

###############################################
################ SHARED SCR ###################
###############################################

loglik.rcm.constant.shared.scr<-function(seroaid.object,par){
	
	loglik<-0
	groups<-unique(seroaid.object$group)
	num.groups<-length(groups)
	age.max<-max(seroaid.object$age)

	log.lambda1<-par[num.groups+1]

	for(i in 1:num.groups){
		
		aux<-which(seroaid.object$group==groups[i])
		
		log.rho1<-par[i]
		
			prob.exp<-sapply(1:age.max,function(x,log.lambda,log.rho)exp(log.lambda)/(exp(log.lambda)+exp(log.rho))*(1-exp(-(exp(log.lambda)+exp(log.rho))*x)),log.lambda=log.lambda1,log.rho=log.rho1)
		
		tabela<-table(factor(seroaid.object$age[aux],levels=1:age.max),factor(seroaid.object$seropos[aux],levels=c(0,1)))	
				
		loglik<-loglik+sum(apply(cbind(tabela[,2],rowSums(tabela),prob.exp),1,function(x)dbinom(x=x[1],size=x[2],prob=x[3],log=T)))
	
				
	} 	
	return(loglik)
}


rcm.constant.estimation.shared.scr<-function(seroaid.object,analysis,n.boot=NULL,conf.level,num.ini){
	
	num.groups<-length(unique(seroaid.object$group))
	groups<-unique(seroaid.object$group)
	age.max<-max(seroaid.object$age)

	cat('Generating initial conditions\n')
	
	output<-rcm.constant.single.group(seroaid.object$age,seroaid.object$seropos,seroaid.object$group,'split',conf.level=conf.level)
	
		
	cat('Performing MLE\n\n')
	
	loglik.max<-(-10^6)
	
	for(n in 1:num.ini){	
		
		estimates.lambda<-log(mean(output$estimates[,2]))
		estimates.rho<-log(output$estimates[,2])+runif(num.groups,-0.1,0.1)
		par<-c(estimates.lambda,estimates.rho)

		cat('Using initial value',n,'out of',num.ini,'\n')	
		sol<-optim(par,loglik.rcm.constant.shared.scr,seroaid.object=seroaid.object,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
		if(sol$value>loglik.max){
			
			sol.best<-sol
			loglik.max<-sol$value
		}
		
	}
	
	sol<-sol.best
	
	fitted.values<-1:age.max	
	estimates<-c()

	lambda1<-exp(sol$par[num.groups+1])

	for(i in 1:num.groups){
		
		rho1<-exp(sol$par[i])

		estimates<-rbind(estimates,c(lambda1,NA,NA,rho1,NA,NA))
		prob.exp<-sapply(1:age.max,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda1,rho=rho1)	

		fitted.values<-cbind(fitted.values,prob.exp)
	}
		
	if(is.null(n.boot)==T){
	
		cat('Calculating CI by inverting LRT\n')
	
		conf.int.lambda<-rcm.constant.estimation.shared.scr.conf.int.lambda(seroaid.object,log(lambda1),sol$value,conf.level=conf.level)

		estimates[,2:3]<-conf.int.lambda$lambda

		cat('Rho\n')
		
		conf.int.rho<-rcm.constant.estimation.shared.scr.conf.int.rho(seroaid.object,log(estimates[,4]),sol$value,conf.level=conf.level)

		estimates[,5:6]<-conf.int.rho$rho

	} else {
			
		cat('Calculating CI by bootstrap\n')
				
		estimates<-rcm.constant.estimation.shared.scr.conf.int.boot(seroaid.object,estimates,n.boot,conf.level=conf.level)

	}
	
	estimates<-data.frame(group=groups,round(estimates,4))

	colnames(estimates)<-c('group','lambda.est','lambda.lower','lambda.upper','rho.est','rho.lower','rho.upper')
	
	colnames(fitted.values)<-c('age',groups)
	
	output<-list(loglik.total=sol$value,estimates=estimates,expected.seroprevalence=fitted.values,df=dim(estimates)[1]+1,model='M1',analysis=analysis)
	
}

####################################################################
###################### Confidence intervals ########################
####################################################################

############################# 1. Lambda ###############################

loglik.rcm.constant.shared.scr.fixed.lambda<-function(seroaid.object,log.lambda,log.rho)loglik.rcm.constant.shared.scr(seroaid.object,par=c(log.rho,log.lambda))		

rcm.constant.shared.scr.estimating.rho.fixed.lambda<-function(seroaid.object,log.lambda){
	
	g<-length(unique(seroaid.object$group))	
	sol<-optim(par=rep(-2,g),loglik.rcm.constant.shared.scr.fixed.lambda,seroaid.object=seroaid.object,log.lambda=log.lambda,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	return(sol)
	
}


rcm.constant.estimation.shared.scr.conf.int.lambda<-function(seroaid.object,log.lambda,loglik.1,conf.level=conf.level){

	cat('Lambda\n')

	f2<-function(log.lambda,seroaid.object,loglik.0)rcm.constant.shared.scr.estimating.rho.fixed.lambda(seroaid.object,log.lambda)$value-loglik.0
	
	num.par<-length(unique(seroaid.object$group))+1
	
	p<-qchisq(conf.level^(1/num.par),1)
	
	#p<-qchisq(0.95,1)
			
	loglik.0<-p/(-2)+loglik.1
		
	##### Calculating lower bound ######

	sol<-rcm.constant.shared.scr.estimating.rho.fixed.lambda(seroaid.object,-10)
		
	if(sol$value>loglik.0){
		
		lower.lambda<-0
				
	} else {
						
		sol.lambda<-uniroot(f2,c(-10,log.lambda),seroaid.object=seroaid.object,loglik.0=loglik.0)	
		lower.lambda<-exp(sol.lambda$root)
	
	}
	
	##### Calculating Upper bound ######
		
	sol<-rcm.constant.shared.scr.estimating.rho.fixed.lambda(seroaid.object,10)
		
	if(sol$value>loglik.0){
		
		upper.lambda<-Inf
		
	} else {

		sol.lambda<-uniroot(f2,c(log.lambda,10),seroaid.object=seroaid.object,loglik.0=loglik.0)	
		upper.lambda<-exp(sol.lambda$root)
	}	
	
	lambda<-c()
	for(i in 1:(num.par-1))lambda<-rbind(lambda,c(lower.lambda,upper.lambda))
	
	return(list(lambda=lambda))
	
}

########################### 2. Rho ###################################

loglik.rcm.constant.shared.scr.fixed.rho<-function(seroaid.object,param,log.rho1,j){
	
	g<-length(unique(seroaid.object$group))

	if(j==1)all.param<-c(log.rho1,param[1:g])
	
	if(j==g)all.param<-c(param[1:(g-1)],log.rho1,param[g])
	
	if(j!=1&j!=g)all.param<-c(param[1:(j-1)],log.rho1,param[j:g])
		
	return(loglik.rcm.constant.shared.scr(seroaid.object,par=all.param))		

}


rcm.constant.shared.scr.estimating.lambda.fixed.rho<-function(seroaid.object,log.rho1,j){
	
	g<-length(unique(seroaid.object$group))	
	sol<-optim(par=rep(-2,g), loglik.rcm.constant.shared.scr.fixed.rho,seroaid.object=seroaid.object, log.rho1=log.rho1,j=j,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	return(sol)
	
}


rcm.constant.estimation.shared.scr.conf.int.rho<-function(seroaid.object,log.rho,loglik.1,conf.level){

	f2<-function(log.rho1,seroaid.object,loglik.0,j)rcm.constant.shared.scr.estimating.lambda.fixed.rho(seroaid.object,log.rho1,j)$value-loglik.0

	g<-length(unique(seroaid.object$group))

	num.par<-g+1
		
	p<-qchisq(conf.level^(1/num.par),1)
	
	#p<-qchisq(0.95,1)
			
	loglik.0<-p/(-2)+loglik.1

	rho<-c()

	for(i in 1:g){
		
		cat('Rho_',i,'\n')

		estimate<-log.rho[i]
	
		##### Calculating lower bound ######

		sol<-rcm.constant.shared.scr.estimating.lambda.fixed.rho(seroaid.object,-10,i)
		
		if(sol$value>loglik.0){
		
			lower.rho<-0
				
		} else {
						
			sol.rho<-uniroot(f2,c(-10,estimate),seroaid.object=seroaid.object,loglik.0=loglik.0,j=i)	
			lower.rho<-exp(sol.rho$root)
	
		}
	
		##### Calculating Upper bound ######
		
		sol<-rcm.constant.shared.scr.estimating.lambda.fixed.rho(seroaid.object,10,i)
		
		if(sol$value>loglik.0){
		
			upper.rho<-Inf
		
		} else {

			sol.rho<-uniroot(f2,c(estimate,10),seroaid.object=seroaid.object,loglik.0=loglik.0,j=i)	
			upper.rho<-exp(sol.rho$root)
		}	
	
		rho<-rbind(rho,c(lower.rho,upper.rho))
	}
	
	return(list(rho=rho))
	
}

###################################################################
################## Confidence intervals (BOOTSTRAP) ###############
###################################################################

rcm.constant.estimation.shared.scr.conf.int.boot<-function(seroaid.object,estimates,n.boot,conf.level){
	
	num.groups<-length(unique(seroaid.object$group))
	
	estimates.bootstrap<-c()
	
	for(i in 1:n.boot){
		
		cat('Analysing bootstrap sample',i,'\n')
		
		seroaid.object.boot<-seroaid.object
			
		groups<-unique(seroaid.object$group)
		
		for(g in 1:num.groups){
		
			aux<-which(seroaid.object.boot$group==groups[g])
			
			lambda<-estimates[g,1]
			rho<-estimates[g,4]
			
			seroaid.object.boot$seropos[aux]<-sapply(seroaid.object.boot$age[aux],function(t,lambda,rho)rbinom(1,1,prob=(lambda/(lambda+rho)*(1-exp(-(lambda+rho)*t)))),lambda=lambda,rho=rho)

			
		}	
		
		#### Generating initial conditions #### 

		par<-log(c(estimates[,4],estimates[1,1]))
			
		sol<-optim(par,loglik.rcm.constant.shared.scr,seroaid.object=seroaid.object.boot,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
		
		estimates.bootstrap<-rbind(estimates.bootstrap,sol$par)
		
	}
	
	for(g in 1:num.groups){
		
		estimates[g,5]<-exp(quantile(estimates.bootstrap[,g],(1-conf.level)/2))
		estimates[g,6]<-exp(quantile(estimates.bootstrap[,g],(1+conf.level)/2))
		estimates[g,2]<-exp(quantile(estimates.bootstrap[,num.groups+1],(1-conf.level)/2))
		estimates[g,3]<-exp(quantile(estimates.bootstrap[,num.groups+1],(1+conf.level)/2))
	
	}
	
	return(estimates)
	
}


############################
######### TWO RCM ##########
############################

rcm.reduction.analysis.overall<-function(age,sero.pos,time.int,time.step,trace){

	times<-seq(time.int[1],time.int[2],time.step)
				
	loglik<-(-1)*(10^6)
		
	tabela<-table(age,sero.pos)
				
	for(i in 1:length(times)){
	
		fit<-optim(par=runif(3,-4,-2),fn= loglik.rcm.reduction.all.param,time.of.change=times[i],age=age,seropos=sero.pos,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
			
		estimates<-fit$par	
												
		loglik.new<-fit$value
	
		if(loglik<loglik.new){
				
			loglik<-loglik.new
			output<-list(loglik.total=loglik,estimates=estimates,time.of.change=times[i],df=3)	
					
		}
		
		if(trace==T)cat('time of change=',times[i],', log.likelihood=',loglik.new,'\n',sep='')
			
	}

	#print(output)	

	return(output)

}

rcm.reduction.analysis.overall.fixed.srr<-function(age,sero.pos,srr,time.int,time.step,trace){

	times<-seq(time.int[1],time.int[2],time.step)
				
	loglik<-(-1)*(10^6)
		
	tabela<-table(age,sero.pos)
				
	for(i in 1:length(times)){
	
		fit<-optim(par=runif(2,-3,-1),fn= loglik.rcm.reduction.lambda.fixed.srr,srr=srr,time.of.change=times[i],age=age,seropos=seropos,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
			
		estimates<-fit$par	
												
		loglik.new<-fit$value
	
		if(loglik<loglik.new){
				
			loglik<-loglik.new
			output<-list(loglik.total=loglik,estimates=estimates,srr=srr,time.of.change=times[i],df=3)	
					
		}
		
		if(trace==T)cat('time of change=',times[i],', log.likelihood=',loglik.new,'\n',sep='')
			
	}
	##print(output)	
	return(output)

}




rcm.reduction.analysis<-function(scm.object,analysis='overall',time.int=c(1,10),time.step=1,srr=F,trace=F,conf.level=0.95){

	age<-scm.object$age
	sero.pos<-scm.object$seropos
			
	results<-c()

	if(analysis=='overall'&srr==F){
		
		output<-rcm.reduction.analysis.overall(age,sero.pos,time.int,time.step,trace)
					
		age.values<-sort(unique(age))	
			
		lambda<-exp(output$estimates[1:2])
				
		rho<-exp(output$estimates[3])
		
		time.of.change<-output$time.of.change
		
		### Confidence interval ####
						
		conf.int<-rcm.reduction.analysis.overall.conf.int(age,sero.pos,log(lambda),log(rho),time.of.change,output$loglik.total,conf.level=conf.level)
		
		output$estimates<-cbind(c(lambda,rho),conf.int)			
		
		fitted.values<-c(1:max(age))		
		
		exp.seroprev<-sapply(1:max(age),seromodel2.time.of.change,change=time.of.change,x=c(lambda,rho))
		
		fitted.values<-cbind(fitted.values,exp.seroprev)
		
		#print(conf.int)				
								
		output<-list(loglik.total=output$loglik.total,estimates=output$estimates,time.of.change=output$time.of.change,df=4,expected.seroprevalence=fitted.values,model='M2',analysis='overall')
				
		return(output)
	
	} else {
		
		if(analysis=='overall'&is.numeric(srr)==T){
		
			output<-rcm.reduction.analysis.overall.fixed.srr(age,sero.pos,srr,time.int,time.step,trace)
					
			age.values<-sort(unique(age))	
			
			lambda<-exp(output$estimates[1:2])
				
			rho<-exp(srr)
		
			time.of.change<-output$time.of.change
		
			### Confidence interval ####
						
			conf.int<-rcm.reduction.analysis.overall.conf.int.fixed.srr(age,sero.pos,log(lambda),log(rho),time.of.change,output$loglik.total,conf.level=conf.level)
							
			output$estimates<-cbind(exp(output$estimates),conf.int)
		
			output$estimates<-rbind(output$estimates,c(rho,NA,NA))
		
			fitted.values<-c(1:max(age))		
		
			exp.seroprev<-sapply(1:max(age),seromodel2.time.of.change,change=time.of.change,x=c(lambda,rho))
		
			fitted.values<-cbind(fitted.values,exp.seroprev)
										
			output<-list(loglik.total=output$loglik.total,estimates=output$estimates,time.of.change=output$time.of.change,df=4,expected.seroprevalence=fitted.values,model='M2',analysis='overall')
				
			return(output)
	
		} else {
			
			cat('ERROR: Group does not exist in the data set\n') 
		
		}
	}
}

seromodel2.time.of.change <- function(age, change,x){

	lambda1  <- x[1] 
	lambda2  <- x[2] 
	rho  	 <- x[3]
	small <- 1E-6
	b0 <- age-change
	theta1 <- lambda1/(lambda1+rho)
	theta2 <- lambda2/(lambda2+rho)
	p_b <- (age > change)*theta1*(1-exp(-(lambda1+rho)*b0)) 
	return((age <= change+small)*(theta2*(1-exp(-(lambda2+rho)*age))) + (age > change+small)*((p_b-theta2)*exp(-(lambda2+rho)*change)+theta2))
}


loglik.scm2.fixed<-function(n1,n0,t,lambda,rho,time.of.change){
	
	p.binom<-as.numeric(sapply(t,seromodel2.time.of.change,change=time.of.change,x=c(lambda,rho)))
	
	data<-cbind(n1,n1+n0)
	data<-cbind(data,p.binom)
	loglik<-as.numeric(apply(data,1,function(x)dbinom.mod(x=x[1],n=x[2],p=x[3])))
	loglik.total<-sum(loglik)

	return(loglik.total)
}


rcm.reduction.analysis.overall.conf.int<-function(age,seropos,lambda,rho,time.of.change,loglik.1,conf.level=conf.level){

	param.aux<-c(lambda,rho)

	num.par<-3
	
	#loglik1<-loglik.rcm.reduction.all.param(age,seropos,param.aux,time.of.change)
	
	#p<-qchisq(conf.level^(1/num.par),1)
	
	p<-qchisq(conf.level,1)
	
	loglik.0<-p/(-2)+loglik.1
	
	output<-c()
	
	for(j in 1:3){
		
		#if(j==1)cat('Calculating CI for lambda_1\n')
		#if(j==2)cat('Calculating CI for lambda_2\n')
		#if(j==3)cat('Calculating CI for rho\n')		
		
		param1<-param.aux[j]
		param<-param.aux[-j]
		
		f2<-function(param1,param,time.of.change,age,seropos,loglik.0,j)rcm.reduction.estimating.one.param(param1,param,time.of.change,age,seropos,j)-loglik.0

				
		##### Calculating lower bound ######
			
		if(param1<(-20)){
		
			lower.bound<-0
		
		} else {
			
			#print(param)
			sol<-rcm.reduction.estimating.one.param(-20,param,time.of.change,age,seropos,j)
				
			if(sol>loglik.0){
		
				lower.bound<-0
				
			} else {
				
				sol<-tryCatch(uniroot(f2,c(-20,param1),age=age,seropos=seropos,loglik.0=loglik.0,j=j,param=param,time.of.change=time.of.change),error=function(e){
					
					sol<-list(root=NA)
					return(sol)
				})	
				lower.bound<-exp(sol$root)
				
			}
		}
		
		#cat('Lower bound=',lower.bound,'\n')
		
		##### Calculating Upper bound ######
		
		sol<-rcm.reduction.estimating.one.param(3,param,time.of.change,age,seropos,j)
		
		if(sol>loglik.0){
		
			upper.bound<-Inf
		
		} else {
		
			sol<-tryCatch(uniroot(f2,c(param1,3),age=age,seropos=seropos,loglik.0=loglik.0,j=j,param=param,time.of.change),error=function(e){
				
				sol<-list(root=NA)
				return(sol)
			})	
			upper.bound<-exp(sol$root)

		}

		#cat('Upper bound=',upper.bound,'\n')
		output<-rbind(output,c(lower.bound,upper.bound))
		
	}
		
	return(output)
	
}

loglik.rcm.reduction.fixed.param<-function(age,seropos,param,param1,time.of.change,j){

	if(j==1)all.param<-c(param1,param)

	if(j==2)all.param<-c(param[1],param1,param[2])
	
	if(j==3)all.param<-c(param,param1)
		
	return(loglik.rcm.reduction.all.param(age,seropos,all.par=all.param,time.of.change=time.of.change))		

}

loglik.rcm.reduction.fixed.param.fixed.srr<-function(age,seropos,srr,param,param1,time.of.change,j){

	if(j==1)all.param<-c(param1,param)

	if(j==2)all.param<-c(param,param1)
		
	output<-loglik.rcm.reduction.all.param(age,seropos,all.param=c(all.param,srr),time.of.change=time.of.change)
			
	return(output)		

}

loglik.rcm.reduction.lambda.fixed.srr<-function(age,seropos,lambda,srr,time.of.change){

	all.param<-c(lambda,srr)	
	return(loglik.rcm.reduction.all.param(age,seropos,all.param=all.param,time.of.change=time.of.change))
	
}



loglik.rcm.reduction.all.param<-function(age,seropos,all.param,time.of.change){
		
	lambda<-exp(all.param[1:2])
	rho<-exp(all.param[3])

	tabela.aux<-table(age,seropos)
		
	n1<-tabela.aux[,2]
	n0<-tabela.aux[,1]
	t<-as.numeric(rownames(tabela.aux))
	
	p.binom<-as.numeric(sapply(t,seromodel2.time.of.change,change=time.of.change,x=c(lambda,rho)))
	data<-cbind(n1,n1+n0)
	data<-cbind(data,p.binom)
	loglik<-as.numeric(apply(data,1,function(x)dbinom(x=x[1],size=x[2],prob=x[3],log=T)))
	loglik.total<-sum(loglik)
	return(loglik.total)
}

rcm.reduction.estimating.one.param<-function(param1,param,time.of.change,age,seropos,j){
	
	fit<-optim(par=runif(2,-3,-1),fn=loglik.rcm.reduction.fixed.param,param1=param1,time.of.change=time.of.change,age=age,seropos=seropos,j=j,control=list(fnscale=-1,maxit=100000,reltol=1e-20))
	return(fit$value)

}

rcm.reduction.analysis.overall.conf.int.fixed.srr<-function(age,seropos,lambda,rho,time.of.change,loglik.1,conf.level=conf.level){

	param.aux<-lambda
	
	num.par<-2
	
	loglik1<-loglik.rcm.reduction.lambda.fixed.srr(age,seropos,param.aux,rho,time.of.change)
				
	output<-c()

	for(j in 1:2){
		
		param1<-param.aux[j]
		param<-param.aux[-j]
		
		f2<-function(param1,param,srr,time.of.change,age,seropos,loglik.0,j)rcm.reduction.estimating.one.param.fixed.srr(param1,param,rho,time.of.change,age,seropos,j)-loglik.0

		#p<-qchisq(conf.level^(1/num.par),1)
	
		p<-qchisq(conf.level,1)
		
		loglik.0<-p/(-2)+loglik1
				
		##### Calculating lower bound ######
			
		sol<-rcm.reduction.estimating.one.param.fixed.srr(-20,param,rho,time.of.change,age,seropos,j)
		
		if(param1<(-20)){
			
			lower.bound<-0
		
		} else {
				
			if(sol>loglik.0){
		
				lower.bound<-0
				
			} else {
						
				sol<-tryCatch(uniroot(f2,c(-20,param1),srr=rho,age=age,seropos=seropos,loglik.0=loglik.0,j=j,param=param,time.of.change=time.of.change),error=function(e){
					
					sol<-list(root=NA)
					return(sol)
					
				})	
				lower.bound<-exp(sol$root)
			
			}
		
		}
		
			
		##### Calculating Upper bound ######
		
		sol<-rcm.reduction.estimating.one.param.fixed.srr(10,param,rho,time.of.change,age,seropos,j)
		
		if(sol>loglik.0){
		
			upper.bound<-Inf
		
		} else {
		
			sol<-tryCatch(uniroot(f2,c(param1,10),srr=rho,age=age,seropos=seropos,loglik.0=loglik.0,j=j,param=param,time.of.change),error=function(e){
				
				sol<-list(root=NA)
				return(sol)
			})	
			upper.bound<-exp(sol$root)

		}

		output<-rbind(output,c(lower.bound,upper.bound))
		
	}
		
	return(output)
	
}

rcm.reduction.estimating.one.param.fixed.srr<-function(param1,param,srr,time.of.change,age,seropos,j){
	
	fit<-optimize(f= loglik.rcm.reduction.fixed.param.fixed.srr,interval=c(-100,10),srr=srr,param1=param1,time.of.change=time.of.change,age=age,seropos=seropos,j=j,maximum=T)
	return(fit$objective)

}

loglik.rcm.reduction.fixed.param.fixed.srr<-function(age,seropos,param,param1,srr,time.of.change,j){
	
	if(j==1)lambda<-c(param1,param)
	
	if(j==2)lambda<-c(param,param1)
				
return(loglik.rcm.reduction.lambda.fixed.srr(age,seropos,lambda=lambda,srr=srr,time.of.change=time.of.change))		

}

lrt.sero.aid<-function(m1,m2){
	
	aux<-which.min(c(m1$df,m2$df))
	
	if(aux==1){
		
		loglik1<-m1$loglik.total
		loglik2<-m2$loglik.total
		df1<-m1$df
		df2<-m2$df
		
	} else {
		
		loglik1<-m2$loglik.total
		loglik2<-m1$loglik.total
		df1<-m2$df
		df2<-m1$df
		
	}
	
	stat<-(-2)*(loglik1-loglik2)
	df<-df2-df1
	p.value<-1-pchisq(stat,df)
	
	cat("Wilks's likelihood ratio test\n")
	cat('Statistic','d.f.','p-value\n',sep='\t')
	cat(stat,df,p.value,'\n',sep='\t')	
	return(p.value)
	
}











