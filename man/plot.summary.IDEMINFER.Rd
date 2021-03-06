% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/im_sace.R
\name{plot.summary.IDEMINFER}
\alias{plot.summary.IDEMINFER}
\title{Plot survivors only and SACE analysis results}
\usage{
\method{plot}{summary.IDEMINFER}(
  x,
  opt = c("pvalue", "effect"),
  by.sace = TRUE,
  delta0 = 0,
  delta1 = 0,
  sace.delta = NULL,
  ...
)
}
\arguments{
\item{x}{A class \code{summary.IDEMSACE} object generated by summary of
\code{IDEMINFER}}

\item{opt}{contour plots of \code{pvalue} or \code{effect}}

\item{by.sace}{Logical value. If True, create a contour plot for given SACE
sensitivity parameter. Otherwise, create a plot for treatment effect for
given imputation sensitivity parameters}

\item{delta0}{Selected treatment arm 0 sensitivity parameters}

\item{delta1}{Selected treatment arm 1 sensitivity parameters}

\item{sace.delta}{Single SACE sensitivity parameter}

\item{...}{Options for \code{plot}}
}
\description{
Generate a plot of survivor only and survivor average causal effect values
}
\details{
The plot function will only generate the contour plot of p-values or
treatment effects on functional outcomes for survivors only analyses.

For SACE analysis, the plot function generates contour plot of line plot
based on the value of \code{by.sace}.
}
\examples{
\dontrun{
rst.abc <- imData(abc, trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
                 y0=NULL, endfml="Y2",
                 trt.label = c("UC+SBT", "SAT+SBT"),
                 cov=c("AGE"), duration=365, bounds=c(0,100));
rst.fit  <- imFitModel(rst.abc);
rst.imp <- imImpAll(rst.fit, deltas=c(-0.25,0,0.25),
                    normal=TRUE, chains = 2, iter = 2000, warmup = 1000);
rst.infer <- imInfer(rst.imp, n.boot = 100, effect.quantiles = c(0.25,0.5,0.75));
rst.survivors <- summary(rst.infer, opt="survivor");
plot(rst.survivors);}

}
