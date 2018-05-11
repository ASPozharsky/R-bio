# 	 PI.R: calculation of probability of identity (PI, ProbId, etc.)
#    Copyright (C) 2018  Alexandr Pozharskiy <aspozharsky@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

#
#	This script implements method described in following papers:
# 		Taberlet P, Luikart G (1999) Non-invasive genetic sampling 
# 			and individual identification. Biological Journal of the Linnaean Society London 68, 41-55;
# 		Waits LP, Luikart G, Taberlet P (2001) Estimating the probability 
# 			of identity among genotypes in natural populations: Cautions and guidelines. Molecular Ecology 10, 249-256.
#
# Important notes:
#		1. Depends on R packages 'adegenet' and 'hierfstat'
#		2. Works with genotypic data in genepop files (*.gen).
#		3. All individuals must be marked as one population. If analysis 
#		of several populations is required, each populations should be written
#		in separate file.
#		4. The script was developed as solution for one particular research
#		case. Application to any other cases may require adjustment and partial 
#		rework, depending on goals.


library('adegenet')
library('hierfstat')
plotPI <- function(datafile, order) {
	ProbIdent <- function(LocusName) 	{ #Eq. 1 in Taberlet, 1999
		b <- LocusName$AlFreq$Freq^4	
		a <- numeric(1)
		for (i in 1:length(LocusName$AlFreq$Freq)) 	{
			for (j in 1:length(LocusName$AlFreq$Freq)) 	{
			if (i!=j) {a <- c(a, (2*LocusName$AlFreq$Freq[i]*LocusName$AlFreq$Freq[j])^2)}}}
			PID <- sum(a)+sum(b)
		return(PID)
						}
	ProbIdent2 <- function(LocusName, SampleSize)	{ #Eq. 2 in Taberlet, 1999
		a2 <- LocusName$AlFreq$Freq^2	
		a3 <- LocusName$AlFreq$Freq^3
		a4 <- LocusName$AlFreq$Freq^4
		PID <- (SampleSize^3*(2*(sum(a2))^2-sum(a4)) - 
		2*SampleSize*(sum(a3)+2*sum(a2)) +
		SampleSize*(9*sum(a2)+2)-6)/((SampleSize-1)*(SampleSize-2)*(SampleSize-3))
		return(PID)
						}

	PI.Sib <- function(LocusName)		{ #Eq. 3
		a <- LocusName$AlFreq$Freq^2
		PID <- 0.25 + 0.5*sum(a) +0.5*sum(a)^2 - 0.25*sum(a^2)
		return(PID)
					}

	PI.Obs <- function(data, Locus)		{ 			# Observed PID. Calculated as sum of probabilities that two randomly chosen individuals
		m <- as.data.frame(table(data[[Locus]]))	# will belong to each genotype
		prob <- sum((m$Freq/nrow(data))^2)		
		return(prob)
					}
	
	PID_reorder <- function(PID, order) {	# Multiplies PIds subsequently according chosen order: decreasing heterozigosity or PId value
		if (order=="hetero") {
			attach(PID)
			Table <- PID[order(Heter, decreasing=TRUE),]
			detach(PID)
			PID_seq <- Table$Prob_Ident
			PID.U_seq <- Table$Prob_Ident.Unb
			PID.S_seq <- Table$Prob_Ident.Sib
			PID.O_seq <- Table$Prob_Ident.Obs	
					} else	{
		if (order=="pid"){
			PID_seq <- sort(PID$Prob_Ident, decreasing=TRUE)
			PID.U_seq <- sort(PID$Prob_Ident.Unb, decreasing=TRUE)
			PID.S_seq <- sort(PID$Prob_Ident.Sib, decreasing=TRUE)
			PID.O_seq <- sort(PID$Prob_Ident.Obs, decreasing=TRUE)
					} else 	print("The ordering parameter must be \"pid\" or \"hetero\" ")
							}
		for (i in 2:nrow(PID))	{
			PID_seq[i] <- PID_seq[i]*PID_seq[i-1]
			PID.U_seq[i] <- PID.U_seq[i]*PID.U_seq[i-1]
			PID.S_seq[i] <- PID.S_seq[i]*PID.S_seq[i-1]
			PID.O_seq[i] <- PID.O_seq[i]*PID.O_seq[i-1]
				}
			Table2 <- data.frame(c(1:nrow(PID)), PID_seq, PID.U_seq, PID.S_seq, PID.O_seq)
			return(Table2)			
		}
		
##### Main section		
	input.gen <- read.genepop(datafile, ncode=3L)
	data.h <- genind2hierfstat(input.gen)
	data.h$pop <- as.numeric(data.h$pop)
	data.h[,1]<-rep(1, times=nrow(data.h))
	data.st <- basic.stats(data.h, diploid=TRUE, digits=4) # compute all required stats
	AllLoci <- list()
	for (i in 1:length(data.st$pop.freq)) {
		AllLoci[[i]] <- assign (paste("Loc",i, sep=""), 
				list (	Names=row.names(data.st$Ho)[i], 
				Het=data.st$Ho[i,1],
				AlFreq=subset(as.data.frame(data.st$pop.freq[[i]]), !is.nan(Freq), select=-Var2)))
						}
				# collect all loci stats in one list
	a <- PI.Obs(data.h, Locus=2)
	# create a table with all PId
	LocNames <- vector(mode="character", length=length(AllLoci))
	Heter <- vector(mode="numeric", length=length(AllLoci))
	Prob_Ident <- vector(mode="numeric", length=length(AllLoci))
	Prob_Ident.Unb <- vector(mode="numeric", length=length(AllLoci))
	Prob_Ident.Sib <- vector(mode="numeric", length=length(AllLoci))
	Prob_Ident.Obs <- vector(mode="numeric", length=length(AllLoci))
	for (LocNum in 1:length(AllLoci)) 	{
		LocNames[LocNum] <- AllLoci[[LocNum]]$Names
		Heter[LocNum] <- AllLoci[[LocNum]]$Het
		Prob_Ident[LocNum] <- ProbIdent(AllLoci[[LocNum]])
		Prob_Ident.Unb[LocNum] <- ProbIdent2(AllLoci[[LocNum]], SampleSize=nrow(input.gen@tab))
		Prob_Ident.Sib[LocNum] <- PI.Sib(AllLoci[[LocNum]])
		Prob_Ident.Obs[LocNum] <- PI.Obs(data.h[,-1], Locus=(LocNum))	
						}
	PID_table <- data.frame(LocNames, Heter, Prob_Ident, Prob_Ident.Unb, Prob_Ident.Sib, Prob_Ident.Obs)
	plotTable <- PID_reorder(PID_table, order=order)
	png(file="PI_plot.png", width=2000, height=2000,res=300)
	plot(x=plotTable[,1],y=plotTable[,2], type="o",xlab="Number of loci", ylab="ProbI", ylim =c(0,0.8), pch=0, lty=2)
	lines(x=plotTable[,1],y=plotTable[,3],type="o", pch=1, lty=2)
	lines(x=plotTable[,1],y=plotTable[,4],type="o", pch=2, lty=2)
	lines(x=plotTable[,1],y=plotTable[,5], type="o", pch=5, lty=1, lwd=1.5)
	legend(x="right",legend=c("ProbI", "ProbI-unbiased", "ProbI-sib", "ProbI-observed"),pch=c(0,1,2,5))
	dev.off()
	write.table(PID_table, file="PI_result.txt",sep="\t")
	write.table(plotTable, file="PI_result.txt",sep="\t", append=TRUE)
}

