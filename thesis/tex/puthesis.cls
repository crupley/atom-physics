%
%  puthesis.cls   2003.04.22  14:49:50   Mark Senn <mds@purdue.edu>
%
%  INDEX: Purdue University thesis document class
%
%  DESCRIPTION:
%
%    This is a LaTeX document class for Purdue University theses.
%
%  USAGE:
%
%    See http://www.ecn.purdue.edu/ECN/Documents/thesis/usage for how to
%    use puthesis.cls.
%
%  HISTORY:
%
%    See http://www.ecn.purdue.edu/ECN/Documents/thesis/history for a
%    record of all changes made to puthesis.cls.
%    2003.02.11 Added
%                   \newcommand{\@@appendixname}{Appendix}
%    2003.02.11 Added
%                   \coversheet{\@@appendixname}
%                 to appendix.
%    2003.02.13 Added
%                   \usepackage{cite}
%    2003.02.21 Changed
%                   \coversheet{\@@appendixname}
%                 to
%                   \coversheet{\uppercase\expandafter{\@@appendixname}}
%               Delete
%                   \clearpage
%                 and
%                   \addtocounter{page}{-1}
%                 from \@@startthebibliography and \@@startvita.
%               Addded
%                   \clearpage
%                 and
%                   \addtocounter{page}{-1}
%                 to \coversheet.
%               Updated puthesis version identification.
%    2003.03.26 Changed
%                   \def\@@makechapterhead#1{\bf\thechapter. #1}
%                 to
%                   \def\@@makechapterhead#1{\thechapter. #1}
%                 for chem option per Ken Yasuda.
%    2003.03.28 Added edci option per Kitty Laird.
%    2003.03.28 Changed
%                   \renewcommand*{\author}[2]%
%                     {
%                       \renewcommand*{\@@TitleAuthor}{#1}%
%                       \renewcommand*{\@@AbstractAuthor}{#2}%
%                     }
%                   to
%                     \def\@@addpoint#1{%
%                       #1%
%                       \def\prev{\relax}%
%                       \@@addthepoint #1\end
%                     }
%                     \def\@@addthepoint#1{%
%                       \ifx #1\end
%                         \if .\prev
%                           \relax
%                         \else
%                           .%
%                         \fi
%                         \let\next=\relax
%                       \else
%                         \def\prev{#1}%
%                         \let\next=\@@addthepoint
%                       \fi
%                       \next
%                     }
%                     \renewcommand*{\author}[2]%
%                       {
%                         \renewcommand*{\@@TitleAuthor}{#1}%
%                         \renewcommand*{\@@AbstractAuthor}{\@@addpoint{#2}}%
%                       }
%                 and changed
%                     \@@AbstractAuthor.
%                   to
%                       % The {  } is to insure the period on the end of
%                       % \@@AbstractAuthor gets treated as the end of a
%                       % sentence.
%                     \@@AbstractAuthor{  }
%                 This was done to make sure a single period always follows the
%               \@@AbstractAuthor, even in the case when \author has a name
%               with a trailing initial, e.g.,
%                     \author{Christopher C. Niessen}{Niessen, Christopher C.}
%               Chris Niessen sent me the above code.  I reformatted it and
%               put it in puthesis.cls.
%    2003.04.01 Added
%                   \otherblankpages
%                   \nootherblankpages
%                   \chapterblankpages
%                   \nochapterblankpages
%                   \coversheets
%                   \nocoversheets
%               commands.
%                 \@@bibstyle was replaced with \@@deptbibstyle and
%               \@@optionbibstyle.  The \@@deptbibstyle (dependent upon
%               department) is used unless \@@optionbibstyle is set,
%               in which case it is used.
%    2003.04.02 Changed
%                   \def\@@makechapterheadspacea{\vspace*{0.5truein}}
%                   \def\@@makechapterheadspaceb{\vspace*{2\baselineskip}}
%                   \def\@@makeschapterheadspacea{\vspace*{0.65625truein}}
%                   \def\@@makeschapterheadspaceb{\vspace*{2\baselineskip}}
%                 to
%                   \def\@@makechapterheadspacea{\vspace*{0.7truein}}
%                   \def\@@makechapterheadspaceb{\vspace*{36pt}}
%                   \def\@@makeschapterheadspacea{\vspace*{0.7truein}}
%                   \def\@@makeschapterheadspaceb{\vspace*{36pt}}
%                 for chem option.
%    2003.04.03 Changed
%                   \setlength{\dblfloatsep}{2\baselineskip}
%                   \setlength{\dbltextfloatsep}{2\baselineskip}
%                   \setlength{\floatsep}{0pt}
%                   \setlength{\intextsep}{2\baselineskip}
%                   \setlength{\textfloatsep}{2\baselineskip}
%                 to
%                   \setlength{\dblfloatsep}{36pt}
%                   \setlength{\dbltextfloatsep}{36pt}
%                   \setlength{\floatsep}{36pt}
%                   \setlength{\intextsep}{36pt}
%                   \setlength{\textfloatsep}{36pt}
%                 for chem option float parameters.
%               Added
%                   \usepackage{rotate}
%                 so \sidewaystable and \sidewaysfigure would work.
%    2003.04.06 Added "lquotation" environment for "cs" option.
%    2003.04.07 Improved list of schools and departments supported.
%    2003.04.17 Dedication is now centered better on page (excluding margins).
%               Changed
%                   \vfil
%                   \@@NotTableOfContents
%                 to
%                   \vfil
%                   \eject
%                   \@@NotTableOfContents
%               Improved list of schools and departments supported.
%    2003.04.18 Improved list of schools and departments supported.
%               The "me" option now does normal size chapter headings per
%                 Taewook Yoo <tyoo@purdue.edu> (tae-wook.yoo.1).  Changed
%                   \def\@@makechapterhead#1{\large\rm\thechapter. \uppercase{#1}}
%                 to
%                   \def\@@makechapterhead#1{\rm\thechapter. \uppercase{#1}}
%    2003.04.22 Single spaced section, subsection, and subsubsection headings
%               for "che" options.
%    2003.04.23 The "eas" option didn't have any corresponding code.
%               Copied "cs" option code to be used for "eas" option.
%               Added "aip" option to use RevTeX's aip.bst bibliography style.
%                 This was mentioned by Physics grad student Brajesh Kumar Rai
%                 <bkumar@physics.purdue.edu> (brajesh.k.rai.1).
%               Changed
%                   \def\@@startthebibliography{}
%                 to
%                   \def\@@startthebibliography{\coversheet{\@@bibname}}
%                 for "ece" option per Jim Pollaro <pollaro@purdue.edu>
%                 (james.r.pollaro.1).
%               Renamed "even" to "next" for first parameter to \@@other.
%                 Note, I don't think this is used anyway.
%               Added a fourth argument to \@@other for the vertical space
%                 after the heading to print.  This is to put three blank
%                 lines after the "LIST OF REFERENCES" heading per Jim Pollaro
%                 <pollaro@purdue.edu> (james.r.pollaro.1).

\NeedsTeXFormat{LaTeX2e}  % why settle for less?
\ProvidesClass{puthesis}[2003/04/23 Purdue thesis class]
\usepackage{cite}
\usepackage{ifthen}
\newcommand{\ifthen}[2]{\ifthenelse{#1}{#2}{}}
\usepackage{endnotes}
\usepackage{pulongtable}
\usepackage{rotating}
\usepackage{lscape}
\makeatletter
\newcount{\@@tcount}
\newlength{\@@tlength}
\newcommand{\@@dept}{unknown}
\newcommand{\set@@dept}[1]{\renewcommand{\@@dept}{#1}}
\DeclareOption{aae}{\set@@dept{aae}}
\DeclareOption{agecon}{\set@@dept{agecon}}
\DeclareOption{ce}{\set@@dept{ce}}
\DeclareOption{che}{\set@@dept{che}}
\DeclareOption{chem}{\set@@dept{chem}}
\DeclareOption{cs}{\set@@dept{cs}}
\DeclareOption{eas}{\set@@dept{eas}}
\DeclareOption{edci}{\set@@dept{edci}}
\DeclareOption{ece}{\set@@dept{ece}}
\DeclareOption{math}{\set@@dept{math}}
\DeclareOption{me}{\set@@dept{me}}
\DeclareOption{mse}{\set@@dept{mse}}
\DeclareOption{physics}{\set@@dept{physics}}
\DeclareOption{stat}{\set@@dept{stat}}
\newcommand{\pendnotes}{}
\DeclareOption{endnotes}{%
  \let\footnote=\endnote
  \renewcommand{\pendnotes}{%
    \newpage
    \begingroup
%% !  \renewcommand{\baselinestretch}{0.9}
      \setlength{\parindent}{0pt}
      \setlength{\parskip}{1.5ex}
      \renewcommand{\enotesize}{\normalsize}
      \theendnotes
    \endgroup
  }
}
\newcommand{\@@type}{thesis}
\newcommand{\set@@type}[1]{\renewcommand{\@@type}{#1}}
\DeclareOption{bypass}{\set@@type{bypass}}
\DeclareOption{prelim}{\set@@type{prelim}}
\DeclareOption{thesis}{\set@@type{thesis}}
\newboolean{@@uglyheadings}
\setboolean{@@uglyheadings}{false}
\DeclareOption{uglyheadings}{\setboolean{@@uglyheadings}{true}}
\newboolean{@@unset}
\newcommand{\@@optionbibstyle}{}
\newcommand{\set@@optionbibstyle}[1]{\renewcommand{\@@optionbibstyle}{#1}}
\DeclareOption{abbrvnat}{\set@@optionbibstyle{abbrvnat}}
\DeclareOption{aip}{\set@@optionbibstyle{aip}}
\DeclareOption{alpha}{\set@@optionbibstyle{alpha}}
\DeclareOption{apalike}{\set@@optionbibstyle{apalike}}
\DeclareOption{chicago}{\set@@optionbibstyle{chicago}}
\DeclareOption{ieee}{\set@@optionbibstyle{ieee}}
\DeclareOption{ieeetr}{\set@@optionbibstyle{ieeetr}}
\DeclareOption{kluwer}{\set@@optionbibstyle{kluwer}}
\DeclareOption{plain}{\set@@optionbibstyle{plain}}
\DeclareOption{plainnat}{\set@@optionbibstyle{plainnat}}
\DeclareOption{unsrt}{\set@@optionbibstyle{unsrt}}
\DeclareOption{unsrtnat}{\set@@optionbibstyle{unsrtnat}}
\newif\if@openbib
\@openbibfalse
%
\newboolean{@@otherblankpages}
\setboolean{@@otherblankpages}{false}
\newcommand{\otherblankpages}{\setboolean{@@otherblankpages}{true}}
\newcommand{\nootherblankpages}{\setboolean{@@nootherblankpages}{false}}
\DeclareOption{otherblankpages}{\otherblankpages}
\DeclareOption{nootherblankpages}{\nootherblankpages}
%
\newboolean{@@chapterblankpages}
\setboolean{@@chapterblankpages}{true}
\newcommand{\chapterblankpages}{\setboolean{@@chapterblankpages}{true}}
\newcommand{\nochapterblankpages}{\setboolean{@@chapterblankpages}{false}}
\DeclareOption{chapterblankpages}{\chapterblankpages}
\DeclareOption{nochapterblankpages}{\nochapterblankpages}
%
\newboolean{@@coversheets}
\setboolean{@@coversheets}{true}
\newcommand{\coversheets}{\setboolean{@@coversheets}{true}}
\newcommand{\nocoversheets}{\setboolean{@@coversheets}{false}}
\DeclareOption{coversheets}{\coversheets}
\DeclareOption{nocoversheets}{\nocoversheets}
%
\DeclareOption{miser}{
  \nootherblankpages
  \nochapterblankpages
  \nocoversheets
}
\DeclareOption{nomiser}{
  \otherblankpages
  \chapterblankpages
  \coversheets
}
%\DeclareOptions*(\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions
\newcommand{\IW}[1]{\immediate\write16{#1}}
{
  \catcode`\+=13
  \def+{\space}
  \ifthen{\equal{\@@dept}{unknown}}
  {{
    \IW{}
    \IW{You must specify an option for your school or department, for example,}
    \IW{``documentclass[aae]{puthesis}''.++The available school and department}
    \IW{options are:}
    \IW{++SCHOOL+++++++++++++++++++++++++++++++DEPARTMENT++++++++++++++++++++++OPTION}
    \IW{++Aeronautics and Astronautics+++++++++++++++++++++++++++++++++++++++++aae}
    \IW{++Agriculture+++++++++++++++++++++++++ Agricultural Economics++++++++++agecon}
    \IW{++Civil Engineering++++++++++++++++++++++++++++++++++++++++++++++++++++ce}
    \IW{++Chemical Engineering+++++++++++++++++++++++++++++++++++++++++++++++++che}
    \IW{++Education++++++++++++++++++++++++++++Curriculum and Instruction++++++edci}
    \IW{++Electrical and Computer Engineering++++++++++++++++++++++++++++++++++ece}
    \IW{++Materials Engineering++++++++++++++++++++++++++++++++++++++++++++++++mse}
    \IW{++Mechanical Engineering+++++++++++++++++++++++++++++++++++++++++++++++me}
    \IW{++Science++++++++++++++++++++++++++++++Chemistry+++++++++++++++++++++++chem}
    \IW{++Science++++++++++++++++++++++++++++++Computer Science++++++++++++++++cs}
    \IW{++Science++++++++++++++++++++++++++++++Earth and Atmospheric Sciences++eas}
    \IW{++Science++++++++++++++++++++++++++++++Mathematics+++++++++++++++++++++math}
    \IW{++Science++++++++++++++++++++++++++++++Physics+++++++++++++++++++++++++physics}
    \IW{++Science++++++++++++++++++++++++++++++Statistics++++++++++++++++++++++stat}
    \IW{}
    \stop
  }}
}
\LoadClass[12pt,twoside]{report}[1996/??/??]
{\obeyspaces
  \gdef\@@Identification{{
    \IW{}
    \IW{                                      ||}
    \IW{                                      ||}
    \IW{                                  see || below}
    \IW{                                      ||}
    \IW{                                   \string\\ || //}
    \IW{                                    \string\\||//}
    \IW{                                     \string\\//}
    \IW{+----------------------------------------------------------------------------+}
    \IW{| This document is using the 2003/04/23 version of puthesis.cls.  It was     |}
    \IW{| probably read in because you used \string\documentclass{puthesis}.                |}
    \IW{|                                                                            |}
    \IW{| See http://www.ecn.purdue.edu/ECN/Documents/tex/unix for instructions on   |}
    \IW{| how to run the the latex in /usr/site/ecn/tetex/bin.                       |}
    \IW{|                                                                            |}
    \IW{| Get old versions of this puthesis.cls file from                            |}
    \IW{|    http://www.ecn.purdue.edu/ECN/Documents/thesis                          |}
    \IW{| if you suspect a recent change broke your document.  I only support the    |}
    \IW{| latest version of this file installed in the system directory.             |}
    \IW{|                                                                            |}
    \IW{| This file may be being read from your directory or a system directory.     |}
    \IW{| To find out which is being used search for ``puthesis.cls'' near the       |}
    \IW{| beginning of this log.  If it says                                         |}
    \IW{|    (puthesis.cls                                                           |}
    \IW{| it's from your directory.  If it says                                      |}
    \IW{|    (/package/teTeX/1.0.7/sparc-sunos5.6/share/texmf/tex/local/puthesis.cls |}
    \IW{| it's from the default system directory.                                    |}
    \IW{|                                                                            |}
    \IW{| ON ECN SUN COMPUTERS, if you want to use the sytem ``puthesis.cls'' type   |}
    \IW{|    mv -i puthesis.cls puthesis.cls.save                                    |}
    \IW{| This will rename your puthesis.cls to puthesis.cls.save, the -i option     |}
    \IW{| will ask you if you want to write over an existing puthesis.cls.save       |}
    \IW{| file.                                                                      |}
    \IW{+----------------------------------------------------------------------------+}
    \IW{                                     //\string\\}
    \IW{                                    //||\string\\}
    \IW{                                   // || \string\\}
    \IW{                                      ||}
    \IW{                                  see || above}
    \IW{                                      ||}
    \IW{                                      ||}
    \IW{}
  }}
}
\ifthen{\equal{aae}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable.  }
    \def\@@makechapterhead#1{\large\bf\thechapter. #1}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{}
    \def\@@startvita{\relax}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{agecon}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\rm\uppercase{Chapter \thechapter. #1}}
    \def\@@makechapterheadspacea{\vspace*{0.6truein}}
    \def\@@makechapterheadspaceb{}
    \def\@@makeschapterheadspacea{\vspace*{0.6truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{ce}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable.  }
    \def\@@makechapterhead#1{\large\rm\thechapter. \uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\normalfont}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\normalfont}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{che}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\thechapter. \uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.65625truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.65625truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\renewcommand{\baselinestretch}{1}}
    \def\@@sectionseries{}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\renewcommand{\baselinestretch}{1}}
    \def\@@subsectionseries{\relax}
    \def\@@subsubsectionbaselinestretch{\renewcommand{\baselinestretch}{1}}
    \def\@@subsubsectionseries{\relax}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{chem}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
      % The vertical spacing between floats appearing either at the
      % top or at the bottom of a page for double column floats in
      % two-column page format.
    \setlength{\dblfloatsep}{36pt}
      % The vertical spacing between floats and text, for both top
      % and bottom floats for double column floats in two-column page
      % format.
    \setlength{\dbltextfloatsep}{36pt}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
      % The vertical spacing between floats appearing either at the
      % top or at the bottom of a page.
    \setlength{\floatsep}{36pt}
    \def\fnum@table{\tablename~\thetable\\}
      % The vertical spacing above and below a float that appears in
      % the middle of a text page with the h placement argument.
    \setlength{\intextsep}{36pt}
    \def\@@makechapterhead#1{\thechapter. #1}
    \def\@@makechapterheadspacea{\vspace*{0.7truein}}
    \def\@@makechapterheadspaceb{\vspace*{36pt}}
    \def\@@makeschapterheadspacea{\vspace*{0.7truein}}
    \def\@@makeschapterheadspaceb{\vspace*{36pt}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
      % The vertical spacing between floats and text, for both top
      % and bottom floats.
    \setlength{\textfloatsep}{36pt}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{cs}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \newenvironment{lquotation}{\begin{quotation}\renewcommand{\baselinestretch}{1}\reset@font}{\end{quotation}}
    \def\@@makechapterhead#1{\thechapter\quad\uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{eas}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \newenvironment{lquotation}{\begin{quotation}\renewcommand{\baselinestretch}{1}\reset@font}{\end{quotation}}
    \def\@@makechapterhead#1{\thechapter\quad\uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{ece}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{- \thepage\ -}\hfil}
    \def\@@makechapterhead#1{
      \ifthenelse{\boolean{@@inappendix}}
        {\large\bf APPENDIX \thechapter\\\uppercase{#1}}
        {\large\bf\thechapter. \uppercase{#1}}
    }
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{- \thepage\ -}\hfil}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\relax}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{edci}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{apalike}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\large\bf\thechapter. #1}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{}
    \def\@@startvita{\relax}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{math}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\large\bf\thechapter. #1}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{}
    \def\@@startvita{\relax}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{me}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable.  }
    \def\@@makechapterhead#1{\rm\thechapter. \uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
    \nochapterblankpages
  }
\ifthen{\equal{mse}{\@@dept}}
  {
    \def\AppendixFigure{\addtocontents{lof}{Appendix Figure\hfil}}
    \def\AppendixTable{\addtocontents{lot}{Appendix Table\hfil}}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vspace*{0.53125truein}}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{
      \rm
      \ifthen{\boolean{@@inappendix}}
        {APPENDIX }
      \thechapter. \uppercase{#1}
    }
    \def\@@makechapterheadspacea{\vspace*{0.65625truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.65625truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\renewcommand*{\baselinestretch}{1.0}}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\renewcommand*{\baselinestretch}{1.0}}
    \def\@@subsectionseries{\normalfont}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\normalfont}
    \def\@@thebibliographyparsep{\parsep 6pt}
  }
\ifthen{\equal{physics}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\large\bf\thechapter. #1}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\ifthen{\equal{stat}{\@@dept}}
  {
    \def\AppendixFigure{\relax}
    \def\AppendixTable{\relax}
    \def\@@bibname{LIST OF REFERENCES}
    \def\@@coversheetspace{\vfill}
    \def\@@deptbibstyle{unsrt}
    \def\@@evenfoot{}
    \def\@@evenhead{\hfil\textrm{\thepage}}
    \def\fnum@table{\tablename~\thetable\\}
    \def\@@makechapterhead#1{\large\rm\thechapter. \uppercase{#1}}
    \def\@@makechapterheadspacea{\vspace*{0.5truein}}
    \def\@@makechapterheadspaceb{\vspace{0.5\baselineskip}}
    \def\@@makeschapterheadspacea{\vspace*{0.5truein}}
    \def\@@makeschapterheadspaceb{\vspace*{0.5\baselineskip}}
    \def\@@oddfoot{}
    \def\@@oddhead{\hfil\textrm{\thepage}}
    \def\@@sectionbaselinestretch{\relax}
    \def\@@sectionseries{\bfseries}
    \def\@@startthebibliography{\coversheet{\@@bibname}}
    \def\@@startvita{\coversheet{\@@vitaname}}
    \def\@@subsectionbaselinestretch{\relax}
    \def\@@subsectionseries{\bfseries}
    \def\@@subsubsectionbaselinestretch{\relax}
    \def\@@subsubsectionseries{\bfseries}
    \def\@@thebibliographyparsep{\relax}
  }
\AtEndDocument{
  \@@Identification
%%  \write\@auxout{\string\newcounter{\string\csname@@@appendix\endcsname}}
%%  \write\@auxout{\string\setcounter{\string\csname@@@appendix\endcsname}{\arabic{@@appendix}}}
}
\newenvironment{cland}
  {\begin{landscape}\hbox\bgroup\hss\vbox\bgroup}
  {\egroup\hss\egroup\end{landscape}}
\newcommand{\coversheet}[1]{
  {
    \ifthen{\boolean{@@coversheets}}
      {
        \ifthenelse{\boolean{@@chapterblankpages}}
          {\clearpage} %% ! {\cleardoublepage}
          {\clearpage}
        \pagestyle{empty}
        \mbox{}
        \@@coversheetspace
        \begin{center}
          #1
        \end{center}
        \vfill
        \newpage
        \addtocounter{page}{-1}
      }
    \@@NotTableOfContents
  }
}
\setlength{\paperheight}{11truein}
\setlength{\paperwidth}{8.5truein}
\renewcommand*{\baselinestretch}{1.5}
\setlength{\evensidemargin}{0.5truein}%!!!
\setlength{\oddsidemargin}{0.5truein}
\setlength{\textheight}{8.75truein}
\setlength{\textwidth}{6truein}
\clubpenalty=10000
\widowpenalty=10000
\displaywidowpenalty=10000
\ifthen{\equal{mse}{\@@dept}}
  {\setlength{\parindent}{0.40625truein}}
%\setlength\parskip{0\p@ \@plus \p@}
\pagenumbering{roman}
\newcommand{\@@NotTableOfContents}{%
  \renewcommand*{\@evenhead}{\@@evenhead}
  \renewcommand*{\@oddhead}{\@@oddhead}
  \renewcommand*{\@evenfoot}{\@@evenfoot}
  \renewcommand*{\@oddfoot}{\@@oddfoot}
}
\@@NotTableOfContents
\setlength{\topmargin}{-0.5truein}
\addtolength{\topmargin}{0.04875truein}
\settoheight{\headheight}{- ivx1 -}
\setlength{\headsep}{0.5truein}
\addtolength{\headsep}{-\headheight}
\addtolength{\headsep}{-0.0625truein}
\newcommand*{\@@TitleAuthor}{\relax}
\newcommand*{\@@AbstractAuthor}{\relax}
  % \@@addpoint and \@@addthepoint are used to insure there is always
  % a single period after the \@@AbstractAuthor, even in the case when
  % the \@@AbstractAuthor ends with a period, e.g.,
  %     \author{Christopher C. Niessen}{Niessen, Christopher C.}
\def\@@addpoint#1{%
  #1%
  \def\prev{\relax}%
  \@@addthepoint #1\end
}
\def\@@addthepoint#1{%
  \ifx #1\end
    \if .\prev
      \relax
    \else
      .%
     \fi
    \let\next=\relax
  \else
    \def\prev{#1}%
    \let\next=\@@addthepoint
  \fi
  \next
}
\renewcommand*{\author}[2]%
  {
    \renewcommand*{\@@TitleAuthor}{#1}%
    \renewcommand*{\@@AbstractAuthor}{\@@addpoint{#2}}%
  }
\newcommand*{\@@TitleDegree}{\relax}
\newcommand*{\@@AbstractDegree}{\relax}
\newcommand*{\@@DegreeMonth}{\relax}
\newcommand*{\@@DegreeYear}{\relax}
\newcommand*{\degree}%
  [4]%
  {%
    \renewcommand*{\@@TitleDegree}{#1}%
    \renewcommand*{\@@AbstractDegree}{#2}%
    \renewcommand*{\@@DegreeMonth}{#3}%
    \renewcommand*{\@@DegreeYear}{#4}%
  }
\newcommand*{\@@MajorProf}{\relax}
\newcommand*{\majorprof}[1]{\renewcommand*{\@@MajorProf}{Major Professor: #1}}
\newcommand*{\majorprofs}[1]{\renewcommand*{\@@MajorProf}{Major Professors: #1}}
\newcommand*{\@@Title}{\relax}
\renewcommand*{\title}[1]{\renewcommand*{\@@Title}{#1}}
\renewcommand*{\maketitle}%
  {{
    \renewcommand*{\baselinestretch}{2} \reset@font
    \begin{titlepage}
    \mbox{}
    \vfil
    \vfil
    \begin{center}
      \uppercase\expandafter{\@@Title}
    \end{center}
    \vfil
    \begin{center}
      \ifthen{\equal{bypass}{\@@type}}{A Master's Bypass Report\\}
      \ifthen{\equal{prelim}{\@@type}}{A Preliminary Report\\}
      \ifthen{\equal{thesis}{\@@type}}{A Thesis\\}
      Submitted to the Faculty\\
      of\\
      Purdue University\\
      by\\
      \@@TitleAuthor
    \end{center}
    \vfil
    \begin{center}
      In Partial Fulfillment of the\\
      Requirements for the Degree\\
      of\\
      \@@TitleDegree
    \end{center}
    \ifthenelse{\equal{mse}{\@@dept}}
      {}
      {\vfil}
    \begin{center}
      \@@DegreeMonth\ \@@DegreeYear
    \end{center}
    \ifthenelse{\equal{mse}{\@@dept}}
      {}
      {\vfil\vfil}
    \end{titlepage}
    \setcounter{page}{2}
  }}
\newenvironment{dedication}%
  {%
    \newpage
    \mbox{}
    \vfil
    \begin{center}%
  }%
  {%
    \end{center}%
    \vfil
    \eject
    \@@NotTableOfContents
  }
\newboolean{@@inother}
\setboolean{@@inother}{false}
  % #1 "next" or "odd": start on next or next odd page?
  % #2 what to print at top of page
  % #3 "y" or "n": put in table of contents?
  % #4 amount of extra space to put after heading at top of page
\newcommand{\@@other}[4]{{%
  \@@NotTableOfContents
  \bgroup
    \setboolean{@@inother}{true}
    \renewcommand{\large}{}%
    \renewcommand{\bf}{}%
    \ifthenelse{\equal{ce}{\@@dept}}
      {\chapter*{\uppercase{#2}}}
      {\chapter*{#2}}
    \ifthen{\equal{y}{#3}}
      {
        \ifthenelse{\equal{ce}{\@@dept}}
          {\addcontentsline{toc}{chapter}{\uppercase{#2}}}
          {\addcontentsline{toc}{chapter}{#2}}
      }
  \egroup
  \vspace{#4}
}}
\newenvironment{acknowledgments}%
  {\@@other{next}{ACKNOWLEDGMENTS}{n}{0pt}}%
  {}
\newenvironment{preface}%
  {\@@other{next}{PREFACE}{n}{0pt}}%
  {}
\renewcommand*{\tableofcontents}{
  \@@other{odd}{TABLE OF CONTENTS}{n}{0pt}
  {\leftskip=0pt \noindent\hbox to\textwidth{\hfil Page}\par}
  {%
%   \renewcommand*{\@oddhead}{%
%     \hfil\textrm{\hfil oddhead a \thepage}%
%     \renewcommand*{\@oddhead}{\hfil oddhead b \thepage}%
%   }
%   \renewcommand*{\@evenhead}{%
%     \hfil\textrm{\hfil evenhead a \thepage}%
%     \renewcommand*{\@evenhead}{\hfil evenhead b \thepage}%
%   }
    \output={
      \let \par \@@par
      \ifnum \outputpenalty<-\@M
        \@specialoutput
      \else
        \@makecol
        \@opcol
        \@startcolumn
        \@whilesw \if@fcolmade \fi
          {%
           \@opcol\@startcolumn}%
      \fi
      \ifnum \outputpenalty>-\@Miv
        \ifdim \@colroom<1.5\baselineskip
          \ifdim \@colroom<\textheight
            \@latex@warning@no@line{Text page \thepage\space
              contains only floats}%
            \@emptycol
          \else
            \global \vsize \@colroom
          \fi
        \else
          \global \vsize \@colroom
        \fi
      \else
        \global \vsize \maxdimen
      \fi
      {\leftskip=0pt \noindent\hbox to\textwidth{\hfil Page}\par}
    }
    \renewcommand{\baselinestretch}{1}\reset@font
    \@starttoc{toc}
  }
}
\ifthen{\equal{che}{\@@dept}\or \equal{ece}{\@@dept}\or \equal{me}{\@@dept}}
  {%
    \def\numberline#1{\hb@xt@\@tempdima{#1\hfil}}
    \renewcommand*\l@figure{\@dottedtocline{1}{0em}{2.3em}}
    \renewcommand*\l@table{\@dottedtocline{1}{0em}{2.3em}}
  }
\def\@dottedtocline#1#2#3#4#5{\ifnum #1>\c@tocdepth \else
  \vskip 0.5\baselineskip plus.2\p@
  {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
    \parindent #2\relax\@afterindenttrue
   \interlinepenalty\@M
   \leavevmode
   \@tempdima #3\relax \advance\leftskip \@tempdima \hbox{}\hskip -\leftskip
    #4\nobreak\leaders\hbox{$\m@th \mkern \@dotsep mu.\mkern \@dotsep
       mu$}\hfill \nobreak
           \hbox to\@pnumwidth{\hfil\reset@font\rm #5}\par}\fi}
\renewcommand*{\listoftables}{
  \@@other{next}{LIST OF TABLES}{y}{0pt}
  {\leftskip=0pt \noindent\hbox to\textwidth{Table\hfil Page}\par}
  {%
    \output={
      \let \par \@@par
      \ifnum \outputpenalty<-\@M
        \@specialoutput
      \else
        \@makecol
        \@opcol
        \@startcolumn
        \@whilesw \if@fcolmade \fi
          {%
           \@opcol\@startcolumn}%
      \fi
      \ifnum \outputpenalty>-\@Miv
        \ifdim \@colroom<1.5\baselineskip
          \ifdim \@colroom<\textheight
            \@latex@warning@no@line{Text page \thepage\space
              contains only floats}%
            \@emptycol
          \else
            \global \vsize \@colroom
          \fi
        \else
          \global \vsize \@colroom
        \fi
      \else
        \global \vsize \maxdimen
      \fi
      {\leftskip=0pt \noindent\hbox to\textwidth{Table\hfil Page}\par}
    }
    \renewcommand{\baselinestretch}{1}\reset@font
    \@starttoc{lot}
  }
}
\renewcommand{\listoffigures}{
  \@@other{next}{LIST OF FIGURES}{y}{0pt}
  {\leftskip=0pt \noindent\hbox to\textwidth{Figure\hfil Page}\par}
  {%
    \output={
      \let \par \@@par
      \ifnum \outputpenalty<-\@M
        \@specialoutput
      \else
        \@makecol
        \@opcol
        \@startcolumn
        \@whilesw \if@fcolmade \fi
          {%
           \@opcol\@startcolumn}%
      \fi
      \ifnum \outputpenalty>-\@Miv
        \ifdim \@colroom<1.5\baselineskip
          \ifdim \@colroom<\textheight
            \@latex@warning@no@line{Text page \thepage\space
              contains only floats}%
            \@emptycol
          \else
            \global \vsize \@colroom
          \fi
        \else
          \global \vsize \@colroom
        \fi
      \else
        \global \vsize \maxdimen
      \fi
      {\leftskip=0pt \noindent\hbox to\textwidth{Figure\hfil Page}\par}
    }
    \renewcommand{\baselinestretch}{1}\reset@font
    \@starttoc{lof}
  }
}
\newcommand{\@@startlist}[1]{
  \@@other{odd}{#1}{y}{0pt}%
  \setlength{\LTleft}{\parindent}%
  \setlength{\LTright}{0truein}%
  %% ! change to \setlength[\@@tlength}{\textwidth - \LTleft - \LTright - 2*\tabcolsep - 1truein} later
  \setlength{\@@tlength}{\textwidth}
  \addtolength{\@@tlength}{-\LTleft}
  \addtolength{\@@tlength}{-\LTright}
  \addtolength{\@@tlength}{-\tabcolsep}
  \addtolength{\@@tlength}{-\tabcolsep}
  \addtolength{\@@tlength}{-1truein}
    %% ! Get real width of first column.
  %% ! \setcounter{\@@tcount}{\c@LT@tables}
  %% ! \addtocounter(\@@tcount){1}
  %% ! \def\crap{\csname LT@\romannumeral\@@tcount\endcsname}
  %% ! \showthe\crap
  %% ! \addtolength{\@@tlength}{-?}
}
\newenvironment{symbols}%
  {%
    \@@startlist{SYMBOLS}
    \begin{longtable}{lp{\@@tlength}}%
  }
  {\end{longtable}}
\newenvironment{abbreviations}
  {
    \@@startlist{ABBREVIATIONS}
    \begin{longtable}{lp{\@@tlength}}%
  }
  {\end{longtable}}
\newenvironment{nomenclature}
  {
    \@@startlist{NOMENCLATURE}
    \begin{longtable}{lp{\@@tlength}}%
  }
  {\end{longtable}}
\renewenvironment{glossary}
  {
    \@@startlist{GLOSSARY}
    \begin{longtable}{lp{\@@tlength}}%
  }
  {\end{longtable}}
\renewenvironment{abstract}%
  {%
    \@@other{next}{ABSTRACT}{y}{0pt}
    {%
      \renewcommand*{\baselinestretch}{1.0}\reset@font
      \vspace*{\baselineskip}
      \vbox{
        \noindent
          % The {  } is to insure the period on the end of
          % \@@AbstractAuthor gets treated as the end of a sentence.
        \@@AbstractAuthor{  }
        \@@AbstractDegree,
        Purdue University,
        \@@DegreeMonth,
        \@@DegreeYear.
        {%
          \renewcommand*{\\}{}%
          \@@Title.
        }
        \@@MajorProf.
      }
      \vspace*{\baselineskip}
      \par
    }
  }%
  {\par}
\renewcommand{\chapter}{%
  \ifthenelse{
        \(      \boolean{@@inother} \and \boolean{@@otherblankpages}   \)
    \or \( \not \boolean{@@inother} \and \boolean{@@chapterblankpages} \)
  }
    {\cleardoublepage}
    {\clearpage}
  \@@NotTableOfContents
  \global\@topnum\z@
  \@afterindentfalse
  \secdef \@chapter \@schapter
}
\newboolean{@@inchapters}
\setboolean{@@inchapters}{false}
\def\@chapter[#1]#2{%
  \ifnum \c@secnumdepth >\m@ne
    \refstepcounter{chapter}%
    \typeout{\@chapapp\space\thechapter.}%
    \setboolean{@@unset}{true}
    \ifthen{\equal{agecon}{\@@dept}}
      {\setboolean{@@unset}{false}\addcontentsline{toc}{chapter}{\protect\uppercase{\@chapapp\ \thechapter. #1}}}
    \ifthen{\equal{ce}{\@@dept}}
      {\setboolean{@@unset}{false}\addcontentsline{toc}{chapter}{\protect\numberline{\thechapter}\uppercase{#1}}}
    \ifthen{\equal{mse}{\@@dept}}
      {\setboolean{@@unset}{false}\addcontentsline{toc}{chapter}{\protect\uppercase{\@chapapp\ \thechapter. #1}}}
    \ifthen{\boolean{@@unset}}
      {\addcontentsline{toc}{chapter}{\protect\numberline{\thechapter}#1}}
  \else
    \ifthenelse{\equal{ce}{\@@dept}}
      {\addcontentsline{toc}{chapter}{\uppercase{#1}}}
      {\addcontentsline{toc}{chapter}{#1}}
  \fi
  \chaptermark{#1}%
% \addtocontents{lof}{\protect\addvspace{10\p@}}%
% \addtocontents{lot}{\protect\addvspace{10\p@}}%
  \ifthenelse{\equal{ce}{\@@dept}}
    {\@makechapterhead{\uppercase{#2}}}
    {\@makechapterhead{#2}}
  \@afterheading
  \ifthen{\not \boolean{@@inchapters}}
    {
      \pagenumbering{arabic}%
      \@@inchapterstrue
    }
}
\def\@makechapterhead#1{%
  \@@makechapterheadspacea
  {\centering
    \@@makechapterhead{#1}
    \endgraf
  }
  \@@makechapterheadspaceb
}
\def\@makeschapterhead#1{%
  \@@makeschapterheadspacea
  {\centering
    \ifthenelse{\equal{ce}{\@@dept}}
      {\large\bf\uppercase{#1}}
      {\large\bf #1}
    \endgraf
  }
  \@@makeschapterheadspaceb
}
\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec
    \leavevmode
  \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim \@tempskipa <\z@
    \@tempskipa -\@tempskipa
  \fi
  \if@nobreak
    \everypar{}%
  \else
    \addpenalty{\@secpenalty}%
    \addvspace{\@tempskipa}%
  \fi
  \@ifstar
    {\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}%
}
\newenvironment{cabstract}{\begin{quote}\textbf{Abstract}\quad}{\end{quote}}
\newtheorem{definition}{Definition}[section]
\newtheorem{observation}{Observation}[section]
\newsavebox{\proofbox}
\sbox{\proofbox}{\rule{7pt}{7pt}}
\newtheorem{Proof}{Proof}
\renewcommand{\theProof}{}
\newenvironment{proof}{\begin{Proof}\rm}{\hfill \usebox{\proofbox} \end{Proof}}
\newtheorem{proposition}{Proposition}[section]
\newtheorem{theorem}{Theorem}[section]
\renewcommand{\section}{%
  \@startsection
    {section}%
    {1}%
    {\z@}%
    {24pt}%
    {12pt}%
    {\@@sectionbaselinestretch\reset@font\@@sectionseries}%
}
\ifthen{\boolean{@@uglyheadings}}
  {
    \def\undertext#1{$\underline{\smash{\hbox{#1}}}$}
    \def\section{
      \vspace{36truept}
      \secdef{\@@section}{\@@ssection}
    }
    \def\@@section[#1]#2{
      \refstepcounter{section}
      \sectionmark{#1}
      \hbox to \textwidth{\hss\undertext{\thesection. #2}\hss}\nopagebreak
      \addcontentsline{toc}{section}{\protect\numberline{\csname thesection\endcsname}#1}
      \nopagebreak
    }
    \def\@@ssection#1{
      \sectionmark{#1}
      \hbox to \textwidth{\hss\undertext{#1}\hss}\nopagebreak
      \nopagebreak
    }
  }
\renewcommand{\subsection}{%
  \@startsection
    {subsection}%
    {2}%
    {\z@}%
    {24pt}% \@plus -1ex \@minus -.2ex}%
    {12pt}% 1.5ex \@plus .2ex}%
    {\@@subsectionbaselinestretch\reset@font\@@subsectionseries}%
}
\ifthen{\boolean{@@uglyheadings}}
  {
    \renewcommand{\subsection}[1]{
      \refstepcounter{subsection}
      \vspace{36truept}
      \subsectionmark{#1}
      \hbox to \textwidth{\hss\thesubsection. #1\hss}\nopagebreak
      \addcontentsline{toc}{subsection}{\protect\numberline{\csname thesubsection\endcsname}#1}
      \nopagebreak
    }
  }
\renewcommand{\subsubsection}{
  \@startsection
    {subsubsection}%
    {3}%
    {\z@}%
    {24pt}% \@plus -1ex \@minus -.2ex}%
    {12pt}% 1.5ex \@plus .2ex}%
    {\@@subsubsectionbaselinestretch\reset@font\@@subsubsectionseries}%
}
\ifthen{\equal{mse}{\@@dept}}
  {
    \def\undertext#1{$\underline{\hbox{#1}}$}
    \newbox{\@@Strut}
    \renewcommand{\subsubsection}[1]{
      \refstepcounter{subsubsection}
      \vspace{6truept}
      \subsubsectionmark{#1}
      \noindent \undertext{\vphantom{gjpqy}\thesubsubsection. #1}\newline\nopagebreak
      \addcontentsline{toc}{subsubsection}{\protect\numberline{\csname thesubsubsection\endcsname}#1}%
      \nopagebreak\indent
    }
  }
\ifthen{\boolean{@@uglyheadings}}
  {
    \renewcommand{\subsubsection}[1]{
      \refstepcounter{subsubsection}
      \vspace{36truept}
      \subsubsectionmark{#1}
      \noindent \undertext{\thesubsubsection. #1}\newline\nopagebreak
      \addcontentsline{toc}{subsubsection}{\protect\numberline{\csname thesubsubsection\endcsname}#1}
      \nopagebreak
    }
  }
\ifthen{\boolean{@@uglyheadings}}
  {
    \setcounter{tocdepth}{4}
    \newcommand*{\l@subsubsubsection}{\@dottedtocline{4}{11.2em}{5.0em}}
    \def\subsubsubsectionmark#1{}
    \newcounter{subsubsubsection}[subsubsection]
    \renewcommand{\thesubsubsubsection}{\thesubsubsection.\arabic{subsubsubsection}}
    \newcommand{\subsubsubsection}[1]{
      \refstepcounter{subsubsubsection}
      \vspace{36truept}
      \subsubsubsectionmark{#1}
      \noindent \thesubsubsubsection. #1\newline\nopagebreak
      \addcontentsline{toc}{subsubsubsection}{\protect\numberline{\csname thesubsubsubsection\endcsname}#1}
      \nopagebreak
    }
  }
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \edef\@svsec{\csname the#1\endcsname\hskip 1em}%
  \fi
  \ifthen{\equal{che}{\@@dept}}
    {%
      \ifnum #2=3
        \refstepcounter{#1}%
        \edef\@svsec{\csname the#1\endcsname\hskip 1em}%
      \fi
    }%
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6\relax
      \@hangfrom{\hskip #3\relax\@svsec}%
      {\interlinepenalty \@M #8\par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth
      \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7%
    }%
  \else
    \def\@svsechd{%
      #6%
      \hskip #3\relax
      \@svsec #8\csname #1mark\endcsname
      {#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth
        \else
          \protect\numberline{\csname the#1\endcsname}
        \fi
        #7%
      }%
    }%
  \fi
  \@xsect{#5}%
}
\renewcommand*{\l@chapter}{\@dottedtocline{0}{0em}{1.4em}}
\ifthen{\equal{ece}{\@@dept}}
  {\renewcommand{\figurename}{Fig.}}
\renewcommand{\fnum@figure}{\figurename~\thefigure.}
\newboolean{@@centercaption}
\renewcommand{\caption}{%
  \setboolean{@@centercaption}{true}
  \refstepcounter\@captype \@dblarg{\@caption\@captype}%
}
\newcommand\bcaption{%
  \setboolean{@@centercaption}{false}
  \refstepcounter\@captype \@dblarg{\@caption\@captype}%
}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1 #2}%
  \ifthenelse{\boolean{@@centercaption}}  % center caption
    {
      \begin{center}
        \renewcommand{\baselinestretch}{1.0}\reset@font
        #1 #2%
      \end{center}
    }
    {
      \setlength{\topsep}{0pt}
      \setlength{\parskip}{0pt}
      \setlength{\partopsep}{0pt}
      \begin{quote}
        \renewcommand{\\}{}
        \renewcommand{\baselinestretch}{1.0}\reset@font
        \ifthenelse{\equal{aae}{\@@dept}}
          {#1. #2\par}
          {
            \ifthenelse{\equal{table}{\@captype}}
              {\hfil\strut #1\hfil\break #2\par}
              {#1 #2\par}
          }
      \end{quote}
    }
  \vskip\belowcaptionskip
}
\setlength\belowcaptionskip{5\p@}
\long\def\@makefntext#1{%
  \baselineskip=12pt
  \noindent
  \@makefnmark #1%
}
\ifthenelse{\equal{}{\@@optionbibstyle}}
  {\newcommand{\@@bibstyle}{\@@deptbibstyle}}
  {\newcommand{\@@bibstyle}{\@@optionbibstyle}}
\ifthen{\equal{abbrvnat}{\@@bibstyle}}
  {
    \usepackage{natbib}
    \renewcommand{\bibname}{{\normalsize\rm LIST OF REFERENCES}}
    \bibliographystyle{abbrvnat}
  }
\ifthen{\equal{aip}{\@@bibstyle}}
  {
    \usepackage{revtex}
    \renewcommand{\bibname}{{\normalsize\rm LIST OF REFERENCES}}
    \bibliographystyle{aip}
  }
\ifthen{\equal{alpha}{\@@bibstyle}}
  {\bibliographystyle{alpha}}
\ifthen{\equal{apalike}{\@@bibstyle}}
  {\bibliographystyle{apalike}}
\ifthen{\equal{chicago}{\@@bibstyle}}
  {
    \usepackage{chicago}
    \bibliographystyle{chicago}
  }
\ifthen{\equal{ieee}{\@@bibstyle}}
  {\bibliographystyle{ieee}}
\ifthen{\equal{ieeetr}{\@@bibstyle}}
  {\bibliographystyle{ieeetr}}
\ifthen{\equal{kluwer}{\@@bibstyle}}
  {
    \usepackage{harvard}
    \bibliographystyle{kluwer}
  }
\ifthen{\equal{plain}{\@@bibstyle}}
  {\bibliographystyle{plain}}
\ifthen{\equal{plainnat}{\@@bibstyle}}
  {
    \usepackage{natbib}
    \renewcommand{\bibname}{{\normalsize\rm LIST OF REFERENCES}}
    \bibliographystyle{plainnat}
  }
\ifthen{\equal{unsrt}{\@@bibstyle}}
  {\bibliographystyle{unsrt}}
\ifthen{\equal{unsrtnat}{\@@bibstyle}}
  {
    \usepackage{natbib}
    \renewcommand{\bibname}{{\normalsize\rm LIST OF REFERENCES}}
    \bibliographystyle{unsrtnat}
  }
\renewenvironment{thebibliography}[1]
  {
    \@@startthebibliography
    \@@other{odd}{\@@bibname}{y}{24pt}
    \list
      {\@biblabel{\arabic{enumiv}}}%
      {%
        \settowidth\labelwidth{\@biblabel{#1}}%
        \leftmargin\labelwidth
        \advance\leftmargin\labelsep
        \if@openbib
          \advance\leftmargin\bibindent
          \itemindent -\bibindent
          \listparindent \itemindent
          \parsep \z@
        \fi
        \@@thebibliographyparsep
        \usecounter{enumiv}%
        \let\p@enumiv\@empty
        \renewcommand{\theenumiv}{\arabic{enumiv}}%
        \baselineskip=12pt
      }%
    \if@openbib
      \renewcommand{\newblock}{\par}
    \else
      \renewcommand{\newblock}{\hskip .11em \@plus.33em \@minus.07em}%
    \fi
    \sloppy\clubpenalty4000\widowpenalty4000%
    \sfcode`\.=\@m
  }
  {}
\newboolean{@@inappendix}
\newcommand{\@@appendixname}{Appendix}
\renewcommand{\appendix}{\par
  \setboolean{@@inappendix}{true}
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \coversheet{\uppercase\expandafter{\@@appendixname}}
  \renewcommand{\@chapapp}{\appendixname}%
  \renewcommand{\thechapter}{\Alph{chapter}}
}
\newcommand{\@@vitaname}{VITA}
\newenvironment{vita}
  {
    \@@startvita
    \@@other{odd}{\@@vitaname}{y}{0pt}%
  }
  {}
\def\@starttoc#1{%
  \begingroup
    \@input{\jobname.#1}%
    \if@filesw \expandafter\newwrite\csname tf@#1\endcsname
      \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax
    \fi
    \global\@nobreakfalse
  \endgroup
}
\newcommand{\num}[1]{\ensuremath{#1}}
\newcommand{\ten}[1]{\ensuremath{{}\cdot 10^{#1}}}
\makeatother
