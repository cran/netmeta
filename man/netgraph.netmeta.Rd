% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netgraph.netmeta.R
\name{netgraph.netmeta}
\alias{netgraph.netmeta}
\title{Network graph}
\usage{
\method{netgraph}{netmeta}(
  x,
  seq = x$seq,
  labels = x$trts,
  cex = 1,
  adj = NULL,
  srt.labels = 0,
  offset = if (!is.null(adj) && all(unique(adj) == 0.5)) 0 else 0.0175,
  scale = 1.1,
  col = gs("col.netgraph"),
  plastic = gs("plastic"),
  thickness = gs("thickness"),
  lwd = 5,
  lwd.min = lwd/2.5,
  lwd.max,
  rescale.thickness,
  dim = "2d",
  rotate = 0,
  highlight = NULL,
  col.highlight = "red2",
  scale.highlight = 1,
  multiarm = gs("multiarm"),
  col.multiarm = NULL,
  alpha.transparency = 0.5,
  points = !missing(cex.points),
  cex.points = 1,
  pch.points = 20,
  col.points = if (length(pch.points) == 1 && pch.points == 21) "black" else "red",
  bg.points = "red",
  points.min,
  points.max,
  rescale.pointsize,
  number.of.studies = gs("number.of.studies"),
  cex.number.of.studies = cex,
  col.number.of.studies = "white",
  bg.number.of.studies = "black",
  pos.number.of.studies = 0.5,
  start.layout = ifelse(dim == "2d", "circle", "eigen"),
  eig1 = 2,
  eig2 = 3,
  eig3 = 4,
  iterate = FALSE,
  tol = 1e-04,
  maxit = 500,
  allfigures = FALSE,
  A.matrix = x$A.matrix,
  N.matrix = sign(A.matrix),
  D.matrix = netdistance(N.matrix),
  xpos = NULL,
  ypos = NULL,
  zpos = NULL,
  figure = TRUE,
  ...
)
}
\arguments{
\item{x}{An object of class \code{netmeta} (mandatory).}

\item{seq}{A character or numerical vector specifying the sequence
of treatments arrangement (anticlockwise if \code{start.layout =
"circle"}).}

\item{labels}{An optional vector with treatment labels.}

\item{cex}{The magnification to be used for treatment labels.}

\item{adj}{One, two, or three values in [0, 1] (or a vector /
matrix with length / number of rows equal to the number of
treatments) specifying the x (and optionally y and z) adjustment
for treatment labels.}

\item{srt.labels}{The character string \code{"orthogonal"} (can be
abbreviated), a single numeric or numerical vector with value(s)
between -180 and 180 specifying the angle to rotate treatment
labels (see Details).}

\item{offset}{Distance between edges (i.e. treatments) in graph and
treatment labels for 2-D plots (value of 0.0175 corresponds to a
difference of 1.75\% of the range on x- and y-axis).}

\item{scale}{Additional space added outside of edges
(i.e. treatments).  Increase this value for larger treatment
labels (value of 1.10 corresponds to an additional space of 10\%
around the network graph).}

\item{col}{A single color (or vector of colors) for lines
connecting treatments (edges) if argument \code{plastic =
FALSE}. Length of the vector must be equal to the number of edges
(see list element 'comparisons' in \code{\link{netmeta}}).}

\item{plastic}{A logical indicating whether the appearance of the
comparisons should be in '3D look' (not to be confused with
argument \code{dim}).}

\item{thickness}{Either a character variable to determine the
method to plot line widths (see Details) or a matrix of the same
dimension and row and column names as argument \code{A.matrix}
with information on line width.}

\item{lwd}{A numeric for scaling the line width of comparisons.}

\item{lwd.min}{Minimum line width in network graph. All connections
with line widths below this values will be set to \code{lwd.min}.}

\item{lwd.max}{Maximum line width in network graph. The connection
with the largest value according to argument \code{thickness}
will be set to this value.}

\item{rescale.thickness}{A logical value or R function to scale the
thickness of lines (see Details).}

\item{dim}{A character string indicating whether a 2- or
3-dimensional plot should be produced, either \code{"2d"} or
\code{"3d"}.}

\item{rotate}{A single numeric with value between -180 and 180
specifying the angle to rotate nodes in a circular network.}

\item{highlight}{A character vector identifying comparisons that
should be marked in the network graph, e.g. \code{highlight =
"treat1:treat2"}.}

\item{col.highlight}{Color(s) to highlight the comparisons given by
\code{highlight}.}

\item{scale.highlight}{Scaling factor(s) for the line width(s) to
highlight the comparisons given by \code{highlight}.}

\item{multiarm}{A logical indicating whether multi-arm studies
should be marked in plot.}

\item{col.multiarm}{Either a function from R package
\bold{colorspace} or \bold{grDevice} to define colors for
multi-arm studies or a character vector with colors to highlight
multi-arm studies.}

\item{alpha.transparency}{The alpha transparency of colors used to
highlight multi-arm studies (0 means transparent and 1 means
opaque).}

\item{points}{A logical indicating whether points should be printed
at nodes (i.e. treatments) of the network graph.}

\item{cex.points, pch.points, col.points, bg.points}{Corresponding
size, type, color, and background color for points. Can be a
vector with length equal to the number of treatments.}

\item{points.min}{Minimum point size. All points with size below
this values will be set to \code{points.min}.}

\item{points.max}{Maximum point size in network graph. The node
with the largest value according to argument \code{cex.points}
will be set to this value.}

\item{rescale.pointsize}{A logical value or R function to scale the
point size (see Details).}

\item{number.of.studies}{A logical indicating whether number of
studies should be added to network graph.}

\item{cex.number.of.studies}{The magnification to be used for
number of studies.}

\item{col.number.of.studies}{Color for number of studies.}

\item{bg.number.of.studies}{Color for shadow around number of
studies.}

\item{pos.number.of.studies}{A single value (or vector of values)
in [0, 1] specifying the position of the number of studies on the
lines connecting treatments (edges). Length of the vector must be
equal to the number of edges.}

\item{start.layout}{A character string indicating which starting
layout is used if \code{iterate = TRUE}. If "circle" (default),
the iteration starts with a circular ordering of the vertices; if
"eigen", eigenvectors of the Laplacian matrix are used,
calculated via generic function \code{\link{eigen}} (spectral
decomposition); if "prcomp", eigenvectors of the Laplacian matrix
are calculated via generic function \code{\link{prcomp}}
(principal component analysis); if "random", a random layout is
used, drawn from a bivariate normal.}

\item{eig1}{A numeric indicating which eigenvector is used as x
coordinate if \code{start = "eigen"} or \code{"prcomp"} and
\code{iterate = TRUE}.  Default is 2, the eigenvector to the
second-smallest eigenvalue of the Laplacian matrix.}

\item{eig2}{A numeric indicating which eigenvector is used as
y-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
\code{iterate = TRUE}.  Default is 3, the eigenvector to the
third-smallest eigenvalue of the Laplacian matrix.}

\item{eig3}{A numeric indicating which eigenvector is used as
z-coordinate if \code{start = "eigen"} or \code{"prcomp"} and
\code{iterate = TRUE}.  Default is 4, the eigenvector to the
fourth-smallest eigenvalue of the Laplacian matrix.}

\item{iterate}{A logical indicating whether the stress majorization
algorithm is carried out for optimization of the layout.}

\item{tol}{A numeric for the tolerance for convergence if
\code{iterate = TRUE}.}

\item{maxit}{An integer defining the maximum number of iteration
steps if \code{iterate = TRUE}.}

\item{allfigures}{A logical indicating whether all iteration steps
are shown if \code{iterate = TRUE}. May slow down computations if
set to \code{TRUE} (especially if \code{plastic = TRUE}).}

\item{A.matrix}{Adjacency matrix (\emph{n}x\emph{n}) characterizing
the structure of the network graph. Row and column names must be
the same set of values as provided by argument \code{seq}.}

\item{N.matrix}{Neighborhood matrix (\emph{n}x\emph{n}) replacing
A.matrix if neighborhood is to be specified differently from node
adjacency in the network graph, for example content-based. Row
and column names must be the same set of values as provided by
argument \code{seq}.}

\item{D.matrix}{Distance matrix (\emph{n}x\emph{n}) replacing
A.matrix and N.matrix if distances should be provided
directly. Row and column names must be the same set of values as
provided by argument \code{seq}.}

\item{xpos}{Vector (\emph{n}) of x coordinates.}

\item{ypos}{Vector (\emph{n}) of y coordinates.}

\item{zpos}{Vector (\emph{n}) of z coordinates.}

\item{figure}{A logical indicating whether network graph should be
shown.}

\item{\dots}{Additional graphical arguments (passed on to
\code{\link{plot.default}}).}
}
\value{
A list containing two data frames with information on nodes and
edges.

\bold{List element 'nodes'}
\item{trts}{Treatment names.}
\item{labels}{Treatment labels.}
\item{seq}{Sequence of treatment labels.}
\item{srt}{String rotation.}
\item{xpos}{Position of treatment / edge on x-axis.}
\item{ypos}{Position of treatment / edge on y-axis.}
\item{zpos}{Position of treatment / edge on z-axis (for 3-D
  plots).}
\item{xpos.labels}{Position of treatment labels on x-axis (for 2-D
  plots).}
\item{ypos.labels}{Position of treatment labels on y-axis (for 2-D
  plots).}
\item{offset.x}{Offset of treatment labels on x-axis (for 2-D
  plots).}
\item{offset.y}{Offset of treatment labels on y-axis (for 2-D
  plots).}
\item{cex}{Point size of treatments / edges.}
\item{col}{Color for points.}
\item{pch}{Point type.}
\item{bg}{Background color for points.}
\item{adj.x}{Adjustment for treatment label on x-axis.}
\item{adj.y}{Adjustment for treatment label on y-axis.}
\item{adj.z}{Adjustment for treatment label on z-axis (for 3-D
  plots).}

\bold{List element 'edges'}
\item{treat1}{Name of first treatment.}
\item{treat2}{Name of second treatment.}
\item{n.stud}{Number of studies directly comparing treatments.}
\item{xpos}{Position of number of studies on x-axis.}
\item{ypos}{Position of number of studies on y-axis.}
\item{adj}{Adjustment of number of studies.}
\item{pos.number.of.studies}{Position of number of studies on
  edge.}
\item{col}{Color for edges.}
}
\description{
This function generates a graph of the evidence network.
}
\details{
This function generates a network graph for an R object created
with \code{\link{netmeta}}. R function \code{\link{plot.default}} is used
to create the network graph.

\subsection{Layout of network graph}{
The network is laid out in the plane, where the nodes in the graph
layout correspond to the treatments and edges display the observed
treatment comparisons. For the default setting, nodes are placed on
a circle.  Other starting layouts are "eigen", "prcomp", and
"random" (Rücker & Schwarzer 2015). If \code{iterate = TRUE}, the
layout is further optimized using the stress majorization
algorithm. This algorithm specifies an 'ideal' distance (e.g., the
graph distance) between two nodes in the plane. In the optimal
layout, these distances are best approximated in the sense of least
squares. Starting from an initial layout, the optimum is
approximated in an iterative process called stress majorization
(Kamada and Kawai 1989, Michailidis and de Leeuw 2001, Hu
2012). The starting layout can be chosen as a circle or coming from
eigenvectors of the Laplacian matrix (corresponding to Hall's
algorithm, Hall 1970), calculated in different ways, or
random. Moreover, it can be chosen whether the iteration steps are
shown (argument \code{allfigures = TRUE}).

An optimized circular presentation which typically has a reduced
(sometimes minimal) number of crossings can be achieved by using
argument \code{seq = "optimal"} in combination with argument
\code{start.layout}. Note, is is not possible of prespecify the
best value for argument \code{start.layout} for any situation as
the result depends on the network structure.
}

\subsection{Definition of line widths}{
Argument \code{thickness} providing the line width of edges
(comparisons) can be a matrix of the same dimension as argument
\code{A.matrix} or any of the following character strings (which
can be abbreviated):
\itemize{
\item Proportional to number of studies comparing two treatments
  (\code{thickness = "number.of.studies"}, default)
\item Proportional to inverse standard error of common effects model
  comparing two treatments (\code{thickness = "se.common"})
\item Proportional to inverse standard error of random effects
  model comparing two treatments (\code{thickness = "se.random"})
\item Weight from common effects model comparing two treatments
  (\code{thickness = "w.common"})
\item Weight from random effects model comparing two treatments
  (\code{thickness = "w.random"})
\item Same line width for all comparisons (\code{thickness =
  "equal"})
}

Only evidence from direct treatment comparisons is considered to
determine the line width if argument \code{thickness} is equal to
any but the last method.

Line widths are determined by argument \code{lwd} if all lines have
the same width. This is possible if either argument \code{thickness
= "equal"}, all pairwise comparisons have the same number of
studies for \code{thickness = "number.of.studies"} or all direct
comparisons are equally precise.

Otherwise, the line width of the thickest line is equal to the
value of argument \code{lwd.max} and all lines with a thickness
below the value of argument \code{lwd.min} are set to this
value. Default for argument \code{lwd.max} is \code{4 * lwd}.

Argument \code{rescale.thickness} can be used to provide a function
to specify the relative line width of edges (comparisons). By
default, the square root function \code{\link[base]{sqrt}} is used
in order to lessen differences in line widths. Argument
\code{rescale.thickness = FALSE} or \code{rescale.thickness = I},
i.e., the identity function \code{\link[base]{I}}, can be used to
not rescale line widths.
}

\subsection{Definition of point sizes}{
Points are printed at nodes (treatments) if argument \code{points =
TRUE} or argument \code{cex.points} is provided.

Point sizes are equal to the value of argument \code{cex.points} if
all points are of equal size. 

Otherwise, the point size of the largest point is equal to the
value of argument \code{points.max} and all points smaller than the
value of argument \code{points.min} are set to this value. The
default for argument \code{points.max} is equal to the largest
value provided in argument \code{cex.points} if this largest value
is below or equal to 25. Otherwise the default is \code{points.max
= 8}.

Argument \code{rescale.pointsize} can be used to provide a function
to specify relative point sizes. Point sizes are not rescaled at
all if they are all equal or the largest \code{cex.points} value is
below or equal to 25. Otherwise, the square root function
\code{\link[base]{sqrt}} is used in order to lessen the differences
in point sizes. Argument \code{rescale.pointsize = FALSE} or
\code{rescale.pointsize = I}, i.e., the identity function
\code{\link[base]{I}}, can be used to not rescale point sizes.
}

\subsection{Other settings}{
Argument \code{srt.labels} can be used to specific the rotation (in
degrees) of the treatment labels. If \code{srt.labels} is equal to
\code{"orthogonal"}, treatment labels are orthogonal to the
circle. If \code{srt.labels} is a single numeric, all labels are
rotated by this degree. If \code{srt.labels} is a numeric vector,
it must be of the same length as the number of treatments and
labels are rotated counter-clockwise starting on the right
side. Finally, if \code{srt.labels} is a named numeric vector, it
must be of the same length as the number of treatments and the
names must be equal to the treatment names (and treatment labels
are rotated according to the specified values).

Further, a couple of graphical parameters can be specified, such as
color and appearance of the edges (treatments) and the nodes
(comparisons), whether special comparisons should be highlighted
and whether multi-arm studies should be indicated as colored
polygons. By default, the \code{\link[colorspace]{sequential_hcl}}
function is used to highlight multi-arm studies.

In order to generate 3-D plots (argument \code{dim = "3d"}), R
package \bold{rgl} is necessary. Note, under macOS the X.Org X
Window System must be available (see
\url{https://www.xquartz.org}).
}
}
\examples{
data(smokingcessation)

# Transform data from arm-based format to contrast-based format
#
pw1 <- pairwise(list(treat1, treat2, treat3),
  event = list(event1, event2, event3), n = list(n1, n2, n3),
  data = smokingcessation, sm = "OR")

# Conduct random effects network meta-analysis
#
net1 <- netmeta(pw1, common = FALSE)

# Network graph with default settings
#
netgraph(net1)

\donttest{
data(Senn2013)

# Generation of an object of class 'netmeta' with reference
# treatment 'plac'
#
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
  data = Senn2013, sm = "MD", reference = "plac")

# Network graph with default settings
#
netgraph(net2)

# Network graph with specified order of the treatments and one
# highlighted comparison
#
trts <- c("plac", "benf", "migl", "acar", "sulf",
  "metf", "rosi", "piog", "sita", "vild")
netgraph(net2, highlight = "rosi:plac", seq = trts)

# Same network graph using argument 'seq' in netmeta function
#
net3 <- netmeta(TE, seTE, treat1, treat2, studlab,
  data = Senn2013, sm = "MD", reference = "plac", seq = trts)
netgraph(net3, highlight = "rosi:plac")

# Network graph optimized, starting from a circle, with multi-arm
# study colored
#
netgraph(net2, start = "circle", iterate = TRUE,
  multiarm = TRUE, col.multiarm = "purple")

# Network graph optimized, starting from a circle, with multi-arm
# study colored and all intermediate iteration steps visible
#
netgraph(net2, start = "circle", iterate = TRUE,
  multiarm = TRUE, col.multiarm = "purple",
  allfigures = TRUE)

# Network graph optimized, starting from Laplacian eigenvectors,
# with multi-arm study colored
#
netgraph(net2, start = "eigen",
  multiarm = TRUE, col.multiarm = "purple")

# Network graph optimized, starting from different Laplacian
# eigenvectors, with multi-arm study colored
#
netgraph(net2, start = "prcomp",
  multiarm = TRUE, col.multiarm = "purple")

# Network graph optimized, starting from random initial layout,
# with multi-arm study colored
#
netgraph(net2, start = "random",
  multiarm = TRUE, col.multiarm = "purple")

# Network graph without plastic look and one highlighted comparison
#
netgraph(net2, plastic = FALSE, highlight = "rosi:plac")

# Network graph with same thickness for all comparisons
#
netgraph(net2, thickness = "equal")

# Network graph with changed labels and specified order of the
# treatments
#
netgraph(net2, seq = c(1, 3, 5, 2, 9, 4, 7, 6, 8, 10),
  labels = LETTERS[1:10])

# Rotate treatment labels (orthogonal to circle)
#
netgraph(net2, srt.labels = "o")

# Network graph in 3-D (opens a new device, where you may rotate and
# zoom the plot using the mouse / the mouse wheel).
# The rgl package must be installed for 3-D plots.
#
if (requireNamespace("rgl", quietly = TRUE))
  netgraph(net2, dim = "3d")
}

}
\references{
Hall KM (1970):
An r-dimensional quadratic placement algorithm.
\emph{Management Science},
\bold{17}, 219--29

Hu Y (2012):
\emph{Combinatorial Scientific Computing}, Chapter Algorithms for
Visualizing Large Networks, pages 525--49.
Chapman and Hall / CRC,  Computational Science.

Kamada T, Kawai S (1989):
An algorithm for drawing general undirected graphs.
\emph{Information Processing Letters},
\bold{31}, 7--15

Krahn U, Binder H, König J (2013):
A graphical tool for locating inconsistency in network meta-analyses.
\emph{BMC Medical Research Methodology},
\bold{13}, 35

Michailidis G, de Leeuw J (2001):
Data visualization through graph drawing.
\emph{Computational Statistics},
\bold{16}, 435--50

Rücker G, Schwarzer G (2016):
Automated drawing of network plots in network meta-analysis.
\emph{Research Synthesis Methods},
\bold{7}, 94--107
}
\seealso{
\code{\link{netmeta}}
}
\author{
Gerta Rücker \email{gerta.ruecker@uniklinik-freiburg.de}, Ulrike
  Krahn \email{ulrike.krahn@bayer.com}, Jochem König
  \email{koenigjo@uni-mainz.de}, Guido Schwarzer
  \email{guido.schwarzer@uniklinik-freiburg.de}
}
\keyword{hplot}
