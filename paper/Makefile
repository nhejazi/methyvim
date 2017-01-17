TITLE = paper
.PHONY : pdf

pdf:
		pandoc $(TITLE).md --latex-engine=xelatex -o $(TITLE).pdf \
			--bibliography $(TITLE).bib
