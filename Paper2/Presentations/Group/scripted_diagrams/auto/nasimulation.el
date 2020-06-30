(TeX-add-style-hook
 "nasimulation"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "tikz"))
 :latex)

