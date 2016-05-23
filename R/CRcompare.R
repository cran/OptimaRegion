#Male.Encap<-read.table("Male_encap.txt", header=TRUE)
#Male.CE<-read.table("Male_CE.txt", header=TRUE)

#### CRcompare.R #####

#' @export
CRcompare<-function(X1,y1,X2,y2,lambda=0.04,responseType='TPS',nosim1and2=200,alpha=0.05,LB1,UB1,triangularRegion1=FALSE, vertex11=NULL, vertex21=NULL, maximization1=TRUE,outputPDFFile1="CR_plot1.pdf",outputOptimafile1="Optima1.txt",LB2,UB2,triangularRegion2=FALSE, vertex12=NULL, vertex22=NULL, maximization2=TRUE,outputPDFFile2="CR_plot2.pdf",outputOptimafile2="Optima2.txt",xlab1and2="Protein eaten (mg)",ylab1and2="Carbohydrates eaten (mg)")
{

  # Computes distribution-free bootstrapped confidence intervals on the mean and median distance between the optima of two different responses. Requires program OptimumRegionTps.R to compute confidence regions on the optima of each response.

  # Usage assuming all default options:
  #out<-CRcompare(X1=X1,y1=y1,X2=X2,y2=y2,LB1=LB1,UB1=UB1,LB2=LB2,UB2=UB2)

  # Arguments

  #       X1--nx2 matrix with the values of the 2 regressors (experimental factors) corresponding to the first response. Note: can have replicates. They will be eliminated by the program and the corresponding y-values averaged
  #       y1--nx1 vector of values for the first response corresponding to X1
  #       X2--nx2 matrix with the values of the 2 regressors (experimental factors) corresponding to the second response. Note: can have replicates. They will be eliminated by the program and the corresponding y-values averaged
  #       y2--nx1 vector of values for the second response corresponding to X2
  #       lambda--penalization parameter (larger implies more smoothing) used to fit the Tps model to both data sets
  #       nosim1and2--number of simulations(default=200) used to find each of the two confidence regions of optima
  #       LB1 and LB2--vector of lower bounds for x (2x1 vector) above which the optimum is sought for the first and second responses, respectively
  #       UB1 and UB2--vector of upper bounds for x (2x1 vector) below which the optimum is sought for the first and second responses, respectively
  #       triangularRegion1 and triangularRegion2--logical: if TRUE it will constrain the maximum points of response 1 (resp., response 2) to lie inside a triangle defined by the coordinates (0,0), and those in "vertex11" (resp: vertex12), and "vertex21" (resp: vertext22), see below (in addition to being constrained to lie inside the region defined by LB1 and UB2 (resp: LB2 and UB2)). NOTE: use TRUE when the treatments form a triangular experimental region in shape. If FALSE, optima will only be constrained to lie inside the rectangular region defined by LB1 and UB2 (resp: LB2 and UB2). Default is FALSE.
  #       vertex11, vertext12---2 times 1 vector with coordinates defining one of the 3 vertices of the triangular region where the first (resp: second) response is being optimized. Must be provided if triangularRegion1 (resp: triangularRegion2) is TRUE (NOTE: vertices numbered clockwise, with vertex0 fixed to (0,0))
  #       vertex21, vertext22--2 times 1 vector with coordinates defining a second  vertex of a triangular region where the first (resp: second) response is being optimized. Must be provided if triangularRegion is TRUE
  #       maximization1,maximization2--logical: if TRUE (default) it maximizes response 1 (resp: response 2) if FALSE it minimizes it
  #       xlab1and2--text label for x axis in both confidence region plots (default: "Protein eaten (mg)")
  #       ylab1and2--text label for y axis in both confidence region plots (default: "Carbohydrates eaten (mg)")
  #       outputPDFFile1 and outPDFFile2--name of the PDF file where the CR plot of the first (resp: second) response is saved (default: "CR_plot.pdf")
  #       outputOptimaFile1 and outputOptimaFile2--name of the text file containing the coordinates of all the simulated optima of the first (resp: second) response


  # Value:

  # dist = vector of distances between pairs of points taken from each set of optima
  # mean,median = mean and median of dist
  # ciMean = 95% confidence interval for the mean of dist using bca bootstrapping
  # ciMEdian = 95% confidence interval for the median of dist using bca bootstrapping

  #last two output values are vectors with 5 columns, containing the signicance level, the next two containing the indices of the order statistics used in the calculations and the final two the calculated endpoints of the CI's.

  # In addition, two PDF files with the CR plots and two text files with the coordinates of each set of optima are created upon completion.

  # Uses: OptimaRegion, spam, boot

  # Written by E. del Castillo, Penn State University, Dept. of Industrial Engineering & Dept. of Statistics
  #         and John Hunt and James Rapkin, University of Exeter, Dept. of Biosciences
  # Version: May 17, 2016
  ##########################################################################################################

  # Load required libraries
  #t<-.libPaths()
  #library("spam", lib.loc=t)
  #library("boot", lib.loc=t)
  #library("OptimaRegion", lib.loc=t)

  if(responseType=='TPS'){ #fit TPS models
    # Run OptRegionTps (from package OptimaRegion) twice
    out1<-OptRegionTps(X=X1,y=y1,nosim=nosim1and2,lambda=lambda,alpha=alpha,LB=LB1,UB=UB1,triangularRegion=triangularRegion1,vertex1=vertex11,vertex2=vertex21,maximization=maximization1,xlab=xlab1and2,ylab=ylab1and2,outputPDFFile=outputPDFFile1,outputOptimaFile=outputOptimafile1)
    out2<-OptRegionTps(X=X2,y=y2,nosim=nosim1and2,lambda=lambda,alpha=alpha,LB=LB2,UB=UB2,triangularRegion=triangularRegion2,vertex1=vertex12,vertex2=vertex22,maximization=maximization2,xlab=xlab1and2,ylab=ylab1and2,outputPDFFile=outputPDFFile2,outputOptimaFile=outputOptimafile2)
  }else if(responseType=='Quad') #fit quadratic polynomails instead
  {
    #Run OptRegionQuad (from package OptimaRegion) twice
    out1<-OptRegionQuad(X=X1,y=y1,nosim=nosim1and2,alpha=alpha,LB=LB1,UB=UB1,triangularRegion=triangularRegion1,vertex1=vertex11,vertex2=vertex21,maximization=maximization1,xlab=xlab1and2,ylab=ylab1and2,outputPDFFile=outputPDFFile1)
    out2<-OptRegionQuad(X=X2,y=y2,nosim=nosim1and2, alpha=alpha,LB=LB2,UB=UB2,triangularRegion=triangularRegion2,vertex1=vertex12,vertex2=vertex22,maximization=maximization2,xlab=xlab1and2,ylab=ylab1and2,outputPDFFile=outputPDFFile2)
  }


  # Form as many pairs of points as in the smallest set and compute their euclidean distances
  Dmatrix<-spam::nearest.dist(x=out1$xin,y=out2$xin,upper=TRUE,delta=1e20)
  dist<-spam::diag(Dmatrix)
  mean<-mean(dist)
  median<-median(dist)
  # Bootstrap the distances between points in each CR
  bmean<-boot::boot(dist,sampleMeans,R=1000)
  bmedian<-boot::boot(dist,sampleMedians,R=1000)
  ciMean<-boot::boot.ci(bmean,type=c("bca"))
  ciMedian<-boot::boot.ci(bmedian,type=c("bca"))
  return(list(dist=dist,mean=mean,median=median,ciMean=ciMean$bca,ciMedian=ciMedian$bca))
}
  sampleMeans<-function(W,i){mean(W[i])}
  sampleMedians<-function(W,i){median(W[i])}
