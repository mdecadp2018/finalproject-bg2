@language python
filename = "report"
import os, platform
import chapter
os.chdir("./../markdown")

if platform.system().lower()=="linux": font = "ukai.ttc"
elif platform.system().lower()=="windows": font = "DFKai-SB"
else: font = "Arial"
fontsize = "12"
margin = "1in"
settingFlag = "--latex-engine=xelatex -V lang=chinese -N --toc --highlight-style kate -V documentclass=report  --filter pandoc-fignos --filter pandoc-tablenos --template=template.tex -V \"CJKmainfont:{0}\" -V fontsize={1}pt -V geometry:margin={2} --bibliography=refer.bib --csl=ieee.csl".format(font, fontsize, margin)
print("---Pandoc---")

# get chapter sequence
chapter_list = ""
for i in chapter.sequence():
    chapter_list += "./paragraph/" + i + " "
    
g.es(chapter_list)

os.system("pandoc cover_and_abstract.md " + chapter_list + " reference.md -o ../pdf/report.pdf {}".format(settingFlag))
g.es("PDF 轉換完畢")
print('-'*12)
