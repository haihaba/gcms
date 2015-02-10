##' Function polygonSelect
##' 
##' Function polygonSelect
##' @export
##' @param coordX
##' @param coordY
##' @return inOut
polygonSelect<-function(coordX,coordY){
  require(sp)
  coordinates <- locator(type = "l")
  coordinates$x[length(coordinates$x + 1)] <- coordinates$x[1]
  coordinates$y[length(coordinates$y + 1)] <- coordinates$y[1]
  newPoly <- Polygon(coordinates)
  inOut <- as.logical(point.in.polygon(coordX, coordY, newPoly@coords[, 1], newPoly@coords[, 2]))
  return(inOut)
}