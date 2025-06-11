FractionCallLog2Ratio<-function(TotalPredBreak,DataSeqTestChrom,TotalVariance,Pos,chr)
{
  MatrixOut<-c()
  for (ll in 1:(length(TotalPredBreak)-1))
  {
    indSeg<-c((TotalPredBreak[ll]+1):TotalPredBreak[ll+1])
    NWindows<-length(indSeg)
    VarianceCopyTest<-TotalVariance*NWindows
    SegValue<-median(DataSeqTestChrom[indSeg])
    SegValueNWindow<-SegValue*NWindows
    if (SegValue>= 0.6)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Multi-Copy Duplication",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    if (SegValue<= -1.1)
    {
      
      MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Double Deletion",round(2*(2^(SegValue))),1,0,SegValue))
      
    }
    
    if (SegValue< 0.6 & SegValue> -1.1)
    {
      Median1CopyTest<- -1
      Median3CopyTest<- 0.58
      
      Segment1CopyTest<-NWindows*Median1CopyTest
      
      Segment3CopyTest<-NWindows*Median3CopyTest
      
      
      IntervalStep1<-seq(1,110,by=1)/100
      
      
      #### Inferring Copy Number 3 ####
      if (SegValue>=0)
      {
        Copy3LikelihoodVec<-c()
        Copy3LikelihoodVecConf2<-c()
        Copy3LikelihoodVecConf1<-c()
        
        for (ii in 1:length(IntervalStep1))
        {
          
          MeanStep1<-((Segment3CopyTest)*IntervalStep1[ii])
          SDStep1<- sqrt(VarianceCopyTest*IntervalStep1[ii])
          Copy3LikelihoodVecConf2[ii]<-dnorm(SegValueNWindow+1.96*SDStep1, mean =MeanStep1, sd = SDStep1,log = FALSE)
          Copy3LikelihoodVecConf1[ii]<-dnorm(SegValueNWindow-1.96*SDStep1, mean =MeanStep1, sd = SDStep1,log = FALSE)
          
          Copy3LikelihoodVec[ii]<-dnorm(SegValueNWindow, mean = MeanStep1,sd = SDStep1 ,log = FALSE)
          
        }
        indLikeliMax<-which.max(Copy3LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy3LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy3LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Duplication",3,FractionEstimate,FractionConf,SegValue))
      }
      
      #### Inferring Copy Number 1 ####
      if (SegValue<0)
      {
        Copy1LikelihoodVec<-c()
        Copy1LikelihoodVecConf2<-c()
        Copy1LikelihoodVecConf1<-c()
        for (ii in 1:length(IntervalStep1))
        {
          
          MeanStep1<-((Segment1CopyTest)*IntervalStep1[ii])
          SDStep1<- sqrt(VarianceCopyTest*IntervalStep1[ii])
          Copy1LikelihoodVecConf2[ii]<-dnorm(SegValueNWindow+1.96*SDStep1, mean =MeanStep1, sd = SDStep1,log = FALSE)
          Copy1LikelihoodVecConf1[ii]<-dnorm(SegValueNWindow-1.96*SDStep1, mean =MeanStep1, sd = SDStep1,log = FALSE)
          
          Copy1LikelihoodVec[ii]<-dnorm(SegValueNWindow, mean = MeanStep1,sd = SDStep1 ,log = FALSE)
        }
        indLikeliMax<-which.max(Copy1LikelihoodVec)
        indLikeliMaxConf1<-which.max(Copy1LikelihoodVecConf1)
        indLikeliMaxConf2<-which.max(Copy1LikelihoodVecConf2)
        
        FractionEstimate<-IntervalStep1[indLikeliMax]
        FractionConf1<-abs(IntervalStep1[indLikeliMaxConf1]-FractionEstimate)
        FractionConf2<-abs(IntervalStep1[indLikeliMaxConf2]-FractionEstimate)
        FractionConf<-max(FractionConf1,FractionConf2)
        MatrixOut<-rbind(MatrixOut,c(chr,Pos[TotalPredBreak[ll]+1],Pos[TotalPredBreak[ll+1]],NWindows,"Deletion",1,FractionEstimate,FractionConf,SegValue))
        
      }
      
    }
    
  }
  MatrixOut
}
