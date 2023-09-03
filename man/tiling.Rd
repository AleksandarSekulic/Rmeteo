\name{tiling}
\alias{tiling}
\title{Tiling \link[raster]{raster} or \link[sp]{Spatial-class} Grid or Pixels object}
\description{Tiling \link[raster]{raster} or \link[sp]{Spatial-class} Grid or Pixels (data frame) object to smaller parts with optional overlap.
}
\usage{
tiling(rast,
       tilesize=500,
       overlapping=50,
       aspoints= FALSE, 
       asfiles=FALSE,
       tilename="tile",
       tiles_folder=paste(getwd(),'tiles',sep='/'),
       parallel.processing=FALSE,
       cpus=6,
       ...)
}
\arguments{
  \item{rast}{\link[terra]{SpatRaster}, SpatialPixels* object, SpatialGrid* object or file path to raster object stored on the disk (can be read via \link[terra]{rast}), for more details see \link[terra]{SpatRaster}. The resolution of the raster should be the same in x and y direction.
}
  \item{tilesize}{integer or vector; tile size in number of cells. Can be a vector of tilesize in x and y direction. Total number of tile cells is \code{tilesize[1]} x \code{tilesize[2]}.}
  \item{overlapping}{integer or vector; overlapping in number of cells. Can be a vector of overlapping in x and y direction.}
  \item{aspoints}{character; Posiible values are \code{sf}, \code{terra} \code{sp}. If specified, tiles are returned in form of points as \link[sf]{sf}, \link[terra]{SpaVector} or \link[sp]{SpatialPointsDataFrame}.}
  \item{asfiles}{boolean; if TRUE tiles are stored on local drive as raster objects.}
  \item{tilename}{character; prefix given to file names}
  \item{tiles_folder}{character; destination folder where tiles will be stored. If doesn't exist, the folder will be created.}
  \item{parallel.processing}{boolean; if TRUE parralel processing is performed via  \link[snowfall]{snowfall-calculation}, sfLapply function.}
  \item{cpus}{integer; number of proccesing units.}
  \item{...}{character; additional arguments for for writing files, see \link[terra]{writeRaster}.}
}
\value{ The list of tiles in \link[terra]{SpatRaster-class} format or in \link[sf]{sf}, \link[terra]{SpaVector} or \link[sp]{SpatialPointsDataFrame} format if \code{aspoints=TRUE}.
}
\author{Milan Kilibarda \email{kili@grf.bg.ac.rs}  }
\seealso{
\code{\link{pred.strk}}
}
\examples{
data(regdata) 
str(regdata@sp) # DEM and TWI data for Serbia at 1 km resolution
dem=regdata@sp['dem']
library(sp)
spplot(dem,  col.regions=bpy.colors() )

str(dem)
# tiling dem in tiles 250x250 with 25 cells overlap
tiles = tiling(dem,tilesize=250,overlapping=25,aspoints=TRUE)
# number of tiles
length(tiles)

image(dem)
points(tiles[[1]] , pch='-' ,col ='green')
points(tiles[[2]], pch='.')

str(tiles[[1]])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

