.lin_impute<-function(t,y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-.lextra(p1=c(t[oindex[1]],y[oindex[1]]),p2=c(t[oindex[2]],y[oindex[2]]),x=t[index.imputed[j]])
    }
  }
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed1<-(index[oindex[i]]:index[oindex[i+1]]) ## index which need to be interpolated + bound
      index.imputed2<-index.imputed1[-c(1,length(index.imputed1))] ## index which need to be interpolated
      len2<-length(index.imputed2)
      for (j in 1:len2)
      {
        yin[index.imputed2[j]]<-.linter(p1=c(t[index.imputed1[1]],y[index.imputed1[1]]),p2=c(t[index.imputed1[length(index.imputed1)]],y[index.imputed1[length(index.imputed1)]]),x=t[index.imputed2[j]])
      }
    }
  } 
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-lextra(p1=c(t[oindex[length(oindex)-1]],y[oindex[length(oindex)]-1]),p2=c(t[oindex[length(oindex)]],y[oindex[length(oindex)]]),x=t[index.imputed[j]])
    }
  }
  yin
}

.linter<-function(p1,p2,x)
{
  x1<-p1[1]; y1<-p1[2]
  x2<-p2[1]; y2<-p2[2]
  m<-abs(x1-x)
  n<-abs(x2-x)
  y<-(m*y2+n*y1)/(m+n)
  y
}

.lextra<-function(p1,p2,x)
{
  x1<-p1[1]; y1<-p1[2]
  x2<-p2[1]; y2<-p2[2]
  m<-abs(x1-x)
  n<-abs(x2-x)
  y<-(m*y2-n*y1)/(m-n)
  y
}


.left_const<-function(y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[1]]
    }
  }
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed<-(index[oindex[i]]:index[oindex[i+1]])
      index.imputed<-index.imputed[-c(1,length(index.imputed))] ## index which need interpolated
      len2<-length(index.imputed)
      for (j in 1:len2)
      {
        yin[index.imputed[j]]<-y[oindex[i]]
      }
    }
  }
  
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[length(oindex)]]
    }
  }
  yin
}


.right_const<-function(y,mindex)
{
  len<- length(y)
  index<-1:len 
  mindex <- sort(mindex)
  oindex <- (1:len)[-mindex]
  yin<-rep(NA,len)
  yin[oindex]<-y[oindex]
  
  if(mindex[1]==1)
  {
    index.imputed<-(1:(oindex[1]-1))
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[1]]
    }
  }
  
  
  for (i in 1:(length(oindex)-1))
  {
    if(oindex[i+1]-oindex[i]>1) 
    {
      index.imputed<-(index[oindex[i]]:index[oindex[i+1]])
      index.imputed<-index.imputed[-c(1,length(index.imputed))] ## index which need interpolated
      len2<-length(index.imputed)
      for (j in 1:len2)
      {
        yin[index.imputed[j]]<-y[oindex[i+1]]
      }
    }
  }
  
  if((len-oindex[length(oindex)])>0)
  {
    index.imputed<-(oindex[length(oindex)]+1):len
    len2<-length(index.imputed)
    for (j in 1:len2)
    {
      yin[index.imputed[j]]<-y[oindex[length(oindex)]]
    }
  }
  yin
}