\name{NEWS}
\title{genepop News}
\encoding{UTF-8}
% cf ?news for the format and usage of this file
%
\section{CHANGES UP TO VERSION 1.0.2}{
\subsection{NEW FEATURES}{
\itemize{
\item The name extension and accessibility of the genepop input file is now checked within the R code.
\item Existence of the output file is checked.
\item More accurate messages when differentation options are run on an input file containing only one 'pop'.
}
}
\subsection{CORRECTED BUGS}{
\itemize{
\item Linkage disequilibrium tests incorrectly indicated 'No data' for tests involving the last locus.
\item The ibd() last line of output in the <>.IBD file mis-reported an uninteresting computation. Both the computation and the line of output have been removed.
\item (R package only) Mantel test always used rank correlation statistic.
}
}
}

