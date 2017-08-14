plot.wfrl <- function(x, minn = 0, col = c("blue","red","green"), vertex.size = 2, vertex.color = c("red", "blue"),
                      edgesThickness = FALSE, zoomThick = 10, ...){
  	if (!inherits(x, "wfrl")) stop("use only with \"wfrl\" object")
  
 if(length(x$lambda1)==1)
  {
    Ac <- as.matrix(x$path[[1]] &x$path[[2]])
	Ah <- as.matrix(x$path[[1]] &!x$path[[2]])
	At <- as.matrix(!x$path[[1]] &x$path[[2]])
	Ao <- Ac|Ah|At

	mdi1 <- max(dim(Ao))
	relations 	<- data.frame(from = which(Ao, arr.ind=TRUE)[,1], to =which(Ao,arr.ind=TRUE)[,2]+ mdi1)
	id        	<- as.numeric(names(table(relations[,1])))[table(relations[,1])>minn]
	wh        	<- unlist(apply(as.matrix(id),1,function(add) which(relations[,1]==add)))
	wh          <- wh[1:min(length(wh), 1000)]

	vert 		<- c(unique(relations[wh,1]),unique(relations[wh,2]))
	cols        <- apply(cbind(relations[wh,1],relations[wh,2]-mdi1),1,function(x) c(Ac[x[1],x[2]],Ah[x[1],x[2]],At[x[1],x[2]]))
	cols 		<- apply(cols,2, function(x){ 
	 if(which(x)==1) return(col[1])
	 if(which(x)==2) return(col[2])
	 if(which(x)==3) return(col[3])
	 })
   
	g 	<- graph.data.frame(relations[wh,], directed=TRUE)#, vertices = vert)


	if(edgesThickness){
  	 wjj <- cbind(relations[wh,1],relations[wh,2]-mdi1)
  	 thick        <- apply(as.matrix(1:dim(wjj)[1]), 1, function(i){
	   if(cols[i] == col[1]) return( (x$regCoef[[1]][wjj[i,1],wjj[i,2]] + x$regCoef[[2]][wjj[i,1],wjj[i,2]])/2)
	   if(cols[i] == col[2]) return( x$regCoef[[1]][wjj[i,1],wjj[i,2]])
	   if(cols[i] == col[3]) return( x$regCoef[[2]][wjj[i,1],wjj[i,2]])  
	 })
	
	 plot(g, edge.arrow.size=.3, vertex.size=vertex.size,vertex.label=NA, 
		 vertex.color= c(rep(vertex.color[1],length(unique(relations[wh,1]))), rep(vertex.color[2],length(unique(relations[wh,2])))),
		 edge.color = cols, edge.width = abs(thick) * zoomThick, ...)
	}
	else
	 plot(g, edge.arrow.size=.3, vertex.size=vertex.size,vertex.label=NA, 
		 vertex.color= c(rep(vertex.color[1],length(unique(relations[wh,1]))), rep(vertex.color[2],length(unique(relations[wh,2])))),
		 edge.color = cols, ...)
  }
  else{
    ask = prod(par("mfcol")) < length(x$lambda1) 
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true) 
    }

    for (i in 1:length(x$lambda1)){

		Ac <- as.matrix(x$path[[i]][[1]] &x$path[[i]][[2]])
	    Ah <- as.matrix(x$path[[i]][[1]] &!x$path[[i]][[2]])
	    At <- as.matrix(!x$path[[i]][[1]] &x$path[[i]][[2]])
	    Ao <- Ac|Ah|At

		relations 	<- data.frame(from = which(Ao, arr.ind=TRUE)[,1], to =which(Ao,arr.ind=TRUE)[,2]+mdi1)
		id        	<- as.numeric(names(table(relations[,1])))[table(relations[,1])>minn]
		wh        	<- unlist(apply(as.matrix(id),1,function(add) which(relations[,1]==add)))
		wh          <- wh[1:min(length(wh), 1000)]

		vert 		<- c(unique(relations[wh,1]),unique(relations[wh,2]))
		cols        <- apply(cbind(relations[wh,1],relations[wh,2]-mdi1),1,function(x) c(Ac[x[1],x[2]],Ah[x[1],x[2]],At[x[1],x[2]]))
		cols 		<- apply(cols,2, function(x){ 
		 if(which(x)==1) return(col[1])
		 if(which(x)==2) return(col[2])
		 if(which(x)==3) return(col[3])
		 })

		g 	<- graph.data.frame(relations[wh,], directed=TRUE, vertices = vert)


		if(edgesThickness){
		 wjj <- cbind(relations[wh,1],relations[wh,2]-mdi1)
		 thick        <- apply(as.matrix(1:dim(wjj)[1]), 1, function(ik){
		   if(cols[i] == col[1]) return( (x$regCoef[[i]][[1]][wjj[ik,1],wjj[ik,2]] + x$regCoef[[i]][[2]][wjj[i,1],wjj[i,2]])/2)
		   if(cols[i] == col[2]) return( x$regCoef[[i]][[1]][wjj[ik,1],wjj[ik,2]])
		   if(cols[i] == col[3]) return( x$regCoef[[i]][[2]][wjj[ik,1],wjj[ik,2]])  
		 })

		 plot(g, edge.arrow.size=.3, vertex.size=vertex.size,vertex.label=NA, 
			 vertex.color= c(rep(vertex.color[1],length(unique(relations[wh,1]))), rep(vertex.color[2],length(unique(relations[wh,2])))),
			 edge.color = cols, edge.width = abs(thick) * zoomThick, ...)
		}
		else
		 plot(g, edge.arrow.size=.3, vertex.size=vertex.size,vertex.label=NA, 
			 vertex.color= c(rep(vertex.color[1],length(unique(relations[wh,1]))), rep(vertex.color[2],length(unique(relations[wh,2])))),
			 edge.color = cols, ...)

  
  }
	
 }

}