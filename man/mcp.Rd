% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcp.R
\name{mcp}
\alias{mcp}
\title{Minimum Convex Hull Polygon}
\usage{
mcp(xy, crs = NULL)
}
\arguments{
\item{xy}{Matrix or Data frame of occurrence coordinates}

\item{crs}{Character of coordinate reference system for minimum convex hull}
}
\description{
Generates a minimum convex polygon (MCP; i.e., convex hull) that
is delineated from occurrence locality coordinates.
This function is a wrapper for `chull()` that makes a SpatialPolygons object.
}
