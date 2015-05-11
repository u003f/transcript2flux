The original method uses genome-scale models of E.coli (Feist et al., 2007: iAF1260.xml) and yeast (Österlund et al, 2013: iTO977.xml). The paper highlights that

"Methods that do not make any assumptions regarding a biological objective (iMAT, Lee–12 and RELATCH*) ... incorrectly predicted a zero growth rate in all cases"

This is because growth has no associated genes / proteins in these models. To alleviate the problem, we extend the models through inclusion of some genes associated with growth. As a first pass, we add the primary DNA polymerases.

For E.coli (iAF1260_DNApoly.xml), this is the DNA polymerase III holoenzyme, with associated genes

b0184 / dnaE
b3701 / dnaN
b0215 / dnaQ
b0470 / dnaX
b0640 / holA
b1099 / holB
b4259 / holC
b4372 / holD
b1842 / holE

As each gene is essential, they are associated with the reaction representing growth via and relationships.

In yeast (iTO977_DNApoly.xml), these are alpha DNA polymerase:primase

YNL102W / pol1
YBL035C / pol12
YIR008C / pri1
YKL045W / pri2

delta DNA polymerase

YDL102W / pol3
YJR006W / pol31
YJR043C / pol32

and epsilon DNA polymerase

YPR175W / dpb2
YBR278W / dpb3
YDR121W / dpb4
YNL262W / pol2
