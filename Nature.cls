\ProvidesClass{nature}[24/2/2004 v1.0]
\typeout{A class for writing preprints for the journal Nature}
\typeout{Written by Peter Czoschke}
\typeout{ }

\LoadClass[10pt,twocolumn]{article} %\LoadClass[12pt]{article}
\RequirePackage{cite}
\RequirePackage{times}
\RequirePackage{fullpage}
\RequirePackage{ifthen}

\renewcommand\@biblabel[1]{#1.} % make labels in bibliobraphy be #
\def\@cite#1#2{$^{\mbox{\tiny #1\if@tempswa, #2\fi}}$} % make citations be superscripts, taken from citesupernumber.sty

% Some style parameters:
\setlength{\parindent}{0.39in}
\setlength{\parskip}{18pt}
\newcommand{\spacing}[1]{\renewcommand{\baselinestretch}{#1}}
%\spacing{2}

% Redefine \maketitle for Nature style:
\def\@maketitle{%
  \newpage\spacing{1}\setlength{\parskip}{12pt}%
    {\huge\noindent\sloppy \textsf{\@title} \par}%
    {\noindent\sloppy \@author}%
}

% Define the affiliations environment, list each institution as an \item:
\newenvironment{affiliations}{%
    \setcounter{enumi}{1}%
    \setlength{\parindent}{0in}%
    \slshape\sloppy%
    \begin{list}{\upshape$^{\arabic{enumi}}$}{%
        \usecounter{enumi}%
        \setlength{\leftmargin}{0in}%
        \setlength{\topsep}{0in}%
        \setlength{\labelsep}{0in}%
        \setlength{\labelwidth}{0in}%
        \setlength{\listparindent}{0in}%
        \setlength{\itemsep}{0ex}%
        \setlength{\parsep}{0in}%
        }
    }{\end{list}\par\vspace{12pt}}

% Redefine the abstract environment to be the first bold paragraph:
\renewenvironment{abstract}{%
    \setlength{\parindent}{0in}%
    \setlength{\parskip}{0in}%
    \bfseries%
    }{\par\vspace{-6pt}}

% Redefine the \section command:
\renewcommand{\section}{\@startsection {section}{1}{0pt}%
    {-6pt}{1pt}%
    {\bfseries}%
    }
\renewcommand{\subsection}{\@startsection {subsection}{2}{0pt}%
    {-0pt}{-0.5em}%
    {\bfseries}*%
    }

% Define the methods environment (use \subsection to separate):
\newenvironment{methods}{%
    \section*{METHODS}%
    \setlength{\parskip}{12pt}%
    }{}

% No heading for References section, but eat up the extra space from \section command:
\renewcommand\refname{\vspace{-48pt}\setlength{\parskip}{0pt}}

% Define the addendum environment for Supplementary Info, Acknowledgements, etc.:
\newenvironment{addendum}{%
    \setlength{\parindent}{0in}%
    \scriptsize%
    \begin{list}{Acknowledgements}{%
        \setlength{\leftmargin}{0in}%
        \setlength{\listparindent}{0in}%
        \setlength{\labelsep}{0em}%
        \setlength{\labelwidth}{0in}%
        \setlength{\itemsep}{3pt}%
        \let\makelabel\addendumlabel}
    }
{\end{list}\normalsize}

\newcommand*{\addendumlabel}[1]{\textbf{#1}\hspace{1em}}
