obj_eu <-
function(X)
	D <-  as.matrix(dist(X,upper=TRUE,diag=TRUE))

