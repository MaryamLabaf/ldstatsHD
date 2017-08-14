#### find triangles from adjacency matrix ####
findTriangles <- function(A1){
      TRI   <- A1%*%A1
      AR    <- which(((TRI>0)*A1)==1,arr.ind=TRUE)

      TRIANG <- NULL
      if(dim(AR)[1]>0){
       if(dim(AR)[1]>1000) print("many triangles are found. Computations can take time \n")
       for(i in 1:dim(AR)[1]){
         print(i)
         da <- as.numeric(names(which(table(c(which(A1[AR[i,1],]!=0),which(A1[AR[i,2],]!=0)))==2)))
         la <- length(da)
         if(la>1)
          TRIANG <- rbind(TRIANG,(cbind(rep(AR[i,1],la),rep(AR[i,2],la),da)))
         else
          TRIANG <- rbind(TRIANG,(c(AR[i,1],AR[i,2],da)))
       }
		TRIANG <- t(apply(TRIANG,1,sort))
        TRIANG <- unique(TRIANG)
       }
      
      return(TRIANG)
}

findAlmTriangles <- function(A1){
      TRI   <- A1%*%A1
      AR    <- which((TRI>0) & A1==0,arr.ind=TRUE)

      TRIANG <- apply(AR,1,function(x)
      {
         da <- as.numeric(names(which(table(c(which(A1[x[1],]!=0),which(A1[x[2],]!=0)))==2)))
         la <- length(da)
         if(la>1)
          cbind(rep(x[1],la),rep(x[2],la),da)
         else
          c(x[1],x[2],da)
      })

      if(class(TRIANG)=="list") TRIANG <- t(apply(do.call(rbind,TRIANG),1,sort))
      else                      TRIANG <- t(apply(t(TRIANG),1,sort))
      
      TRIANG <- unique(TRIANG)
      pp <- apply(TRIANG,1,function(x) any(table(x)==2))
      TRIANG <-TRIANG[!pp,]
      return(TRIANG)
}

