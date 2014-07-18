normv2 <-
function(s){
	
	if(nrow(s) == 1){
		sn  <-  s/sqrt(sum(s^2))
		
	}else
		sn  <-  apply(s,1,function(s) s/sqrt(sum(s^2)))
}

