plot.pcorSim <- function(x, vertex.size = 3,vertex.label = NA, hubsCol = TRUE, ...){

  if (!inherits(x, "pcorSim")) stop("use only with \"pcorSim\" object")
   PATH 		<- x$path
   wh          	<- which(apply(PATH,1,sum)>0)
   PATH        <- PATH[wh,wh]
   gr          <- suppressMessages(graph.adjacency(PATH,mode="undirected"))
   if(hubsCol){
    COL         	<- rep(1,length(wh))
    if(length(x$hubs)>0)
    {
      whHubs <- apply(as.matrix(x$hubs),1,function(j) which(j==wh))
      COL[whHubs] <-2
    }
     suppressMessages(plot(gr,vertex.size = vertex.size, vertex.label = vertex.label, vertex.color= COL, ...))
   }
   else{
     suppressMessages(plot(gr,vertex.size = vertex.size, vertex.label = vertex.label, ...))
   }
   
}
