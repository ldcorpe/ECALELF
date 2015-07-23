#!/bin/bash
OUTDIR=comparison

usage(){
    echo "`basename $0`  tablefiles"
    echo "  -e: no error"
    echo "  --outdir"
}

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o he -l help,column:,outdir:,plusMC,noerror -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

set -- $options

while [ $# -gt 0 ]
do
    case $1 in
	-h|--help) usage; exit 0;;
	--noerror) NOERROR=y;;
	--outdir) OUTDIR=$2; shift;;
	--plusMC) plusMC=y;;
	-e|--noerror) error="-e";;
	(--) shift; break;;
	(-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
	(*) break;;
    esac
    shift
done
var=0
mkdir -p $OUTDIR
files=""
filesA=""
filesB=""
rm $OUTDIR/slide*
for file in $@
do
files="$files $file"
if [ $var -lt 5 ];then
echo "add to A"
filesA="$filesA $file"
else
echo "add to B"
filesB="$filesB $file"
fi
((var++))
done

echo "done adding"
it=0;
array=( 2 3 5 6 8 10 11 15 )
for it in "${array[@]}"
do

# !ds = dollar sign $ (sub later, sed gets confused with escapable chars)
# !bs = backslash sign \ (sub later)
# !us = backslash sign _ (sub later)
case $it in 
    2)
	cName="nEvents"
	;;
    3|4)
	cName="!ds !bsDelta m!ds"
	;;
    5)
	cName="!ds !bsDelta P !ds"
	;;
    6)
	cName="!ds !bssigma!us{CB}!ds"
	;;
    8|9)
	cName="!ds!bsresol !ds"
	;;
    10)
	cName="add. smear."
	;;
    11)
	cName="!ds !bschi^2 !ds."
	;;
    15)
	cName="!ds !bssigma!us{eff} !ds"
	;;
    *)
	exit 1
	;;
esac
echo "0var is $var"

if [ $var -lt 5 ];then
echo "var is $var"
tmpFile=$OUTDIR/slide$it.tex
#./script/compareColumns.sh --column $it $files
./script/compareColumns.sh $error --column $it $files
#echo " ./script/compareColumns.sh -e --column $it $files"
#cat tmp/file.tex
cp tex/compareSlides.tex $tmpFile 
sed -i "/_TABLE_/ r tmp/file.tex" $tmpFile
sed -i '/_TABLE_/ d' $tmpFile
sed -i "s|0poiname|$cName|g " $tmpFile
sed -i "s|!us|_|g " $tmpFile
sed -i "s|!bs|\\\|g " $tmpFile
sed -i "s|!ds|$|g " $tmpFile
else
tmpFileB=$OUTDIR/slide$it.B.tex
tmpFileA=$OUTDIR/slide$it.A.tex
#./script/compareColumns.sh --column $it $files
echo "files A $filesA"
echo "files B $filesB"
./script/compareColumns.sh $error --column $it $filesB
cp tex/compareSlides.tex $tmpFileB
sed -i "/_TABLE_/ r tmp/file.tex" $tmpFileB

sed -i '/_TABLE_/ d' $tmpFileB
sed -i "s|0poiname|$cName|g " $tmpFileB
sed -i "s|!us|_|g " $tmpFileB
sed -i "s|!bs|\\\|g " $tmpFileB
sed -i "s|!ds|$|g " $tmpFileB

./script/compareColumns.sh $error --column $it $filesA
cp tex/compareSlides.tex $tmpFileA 
sed -i "/_TABLE_/ r tmp/file.tex" $tmpFileA

sed -i '/_TABLE_/ d' $tmpFileA
sed -i "s|0poiname|$cName|g " $tmpFileA
sed -i "s|!us|_|g " $tmpFileA
sed -i "s|!bs|\\\|g " $tmpFileA
sed -i "s|!ds|$|g " $tmpFileA


fi

done


cat > $OUTDIR/slide0.tex <<EOF
\begin{comment}

\end{comment}

\makeatletter
\@ifundefined{dataSample}{
  \newcommand{\dataSample}{$dataSample}
  \newcommand{\mcSample}{$mcSample}
  \newcommand{\imgDir}{$dirData/$selection/$invMass_var/img}
  \newcommand{\PeriodDivisionTickzOld}{}
  \newcommand{\PeriodDivisionTickzNew}{}
  \newcommand{\PeriodDivisionTickz}{}
 \newcommand{\period}{$PERIOD}
}{
  %\renewcommand{\dataSample}{$dataSample}
  \renewcommand{\mcSample}{$mcSample}
  \renewcommand{\imgDir}{$dirData/$selection/$invMass_var/img}
  \renewcommand{\period}{$PERIOD}
}
\@ifundefined{invMassVarName}{
  \newcommand{\invMassVarName}{$invMassVarName}
}{
  \renewcommand{\invMassVarName}{$invMassVarName}
}
\@ifundefined{resol}{
  \newcommand{\resol}{\ensuremath{\frac{\sigma_{CB}}{peak_{CB}}}}
}{}

\makeatother

\usebackgroundtemplate{
	\includegraphics[width=\paperwidth,height=\paperheight]{logos/Blasenkammer_by_HenryGale_white}%
}

\section{\invMassVarName}
\frame{\centering\scalebox{2}{\textbf{\textsf{\rotatebox{35}{Rereco comparison tables}}}}}

\usebackgroundtemplate{
	\includegraphics[width=\paperwidth,height=\paperheight]{}%
}



EOF

#cat tex/template.tex $OUTDIR/slide2.tex $OUTDIR/slide3.tex  $OUTDIR/slide5.tex  $OUTDIR/slide6.tex  $OUTDIR/slide8.tex  $OUTDIR/slide10.tex  $OUTDIR/slide11.tex  $OUTDIR/slide15.tex template_end.tex > $OUTDIR/$OUTDIR.tex
cat tex/template.tex $OUTDIR/slide*.tex tex/template_end.tex > $OUTDIR/$OUTDIR.tex


pdflatex $OUTDIR/$OUTDIR.tex
cp $OUTDIR.pdf slides/.
