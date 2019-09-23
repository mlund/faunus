---
title: \emph{The Faunus User Manual}
papersize: a4
fontsize: 10pt
geometry: margin=3cm
toccolor: Maroon
linkcolor: Maroon
citecolor: Maroon
urlcolor: Maroon
# https://tex.stackexchange.com/questions/3214/error-when-using-t1-fontenc-and-urw-garamond-from-mathdesign
# fontfamily: mathdesign
# fontfamilyoptions: urw-garamond
mainfont: EB Garamond       # https://github.com/octaviopardo/EBGaramond12
mathfont: Garamond-Math     # https://github.com/YuanshengZhao/Garamond-Math
# sansfont: Gill Sans Nova
monofont: Fira Code         # https://github.com/tonsky/FiraCode
monofontoptions:
- StylisticSet=4 # sans-serif like g
- StylisticSet=5 # traditional shape of @
# monofont: Inconsolata
# monofont: Source Code Pro
linestretch: 1.15
links-as-notes: true
header-includes: |
    \usepackage{fancyhdr}
    \pagestyle{fancy}
    \fancyhead[CO,CE]{Faunus}
    \fancyfoot[CO,CE]{}
    \fancyfoot[LE,RO]{\thepage}
#    \usepackage[margins=raggedright]{floatrow}

---

